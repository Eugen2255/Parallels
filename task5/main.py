import argparse
import time
import queue
import threading
import gc
import cv2
import statistics
import os
import numpy as np
from typing import Optional, Dict, List
from dataclasses import dataclass
from ultralytics import YOLO


@dataclass
class FrameData:
    frame_id: int
    image: np.ndarray
    result: Optional[np.ndarray] = None
    processed: bool = False

class ResourceRAII:
    def __init__(self, resource_type: str, **kwargs):
        self.resource_type = resource_type
        self.resource = None
        self.kwargs = kwargs
        self._acquire()
    
    def _acquire(self):
        if self.resource_type == "capture":
            self.resource = cv2.VideoCapture(self.kwargs.get("source", 0))
            if not self.resource.isOpened():
                raise RuntimeError(f"Не удалось открыть {self.kwargs.get('source', 'камеру')}")
            self.resource.set(cv2.CAP_PROP_FRAME_WIDTH, 640)
            self.resource.set(cv2.CAP_PROP_FRAME_HEIGHT, 480)
        elif self.resource_type == "writer":
            fourcc = cv2.VideoWriter_fourcc(*'mp4v')
            self.resource = cv2.VideoWriter(
                self.kwargs["filename"],
                fourcc,
                self.kwargs.get("fps", 25.0),
                (self.kwargs.get("width", 640), self.kwargs.get("height", 480))
            )
    
    def __enter__(self):
        return self.resource
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._release()
        return False
    
    def __del__(self):
        self._release()
    
    def _release(self):
        """Освобождение ресурса"""
        if self.resource is not None:
            self.resource.release()
            self.resource = None
    
    def get(self):
        return self.resource


class WorkerThread(threading.Thread):
    def __init__(self, input_queue: queue.Queue, output_queue: queue.Queue, 
                 model_path: str, stop_event: threading.Event):
        super().__init__(daemon=True)
        self.input_queue = input_queue
        self.output_queue = output_queue
        self.model_path = model_path
        self.stop_event = stop_event
        # Модель создаётся внутри потока
        self.model: Optional[YOLO] = None
    
    def run(self):
        self.model = YOLO(self.model_path)
        
        while not self.stop_event.is_set():
            try:
                frame_data: FrameData = self.input_queue.get(timeout=0.1)
            except queue.Empty:
                continue
            
            if frame_data is None:
                self.input_queue.task_done()
                break
            
            try:
                results = self.model.predict(frame_data.image, verbose=False)
                
                if results and results[0].keypoints is not None:
                    frame_data.result = results[0].plot()
                else:
                    frame_data.result = frame_data.image.copy()
                
                frame_data.processed = True
                
            except Exception as e:
                print(f"Ошибка обработки кадра {frame_data.frame_id}: {e}")
                frame_data.result = frame_data.image.copy()
            
            try:
                self.output_queue.put(frame_data, timeout=0.1)
            except queue.Full:
                pass 
            
            self.input_queue.task_done()

class PoseInferencePipeline:
    def __init__(self, video_path: str, output_path: str, 
                 mode: str = "single", num_workers: int = 4,
                 model_path: str = "yolov8s-pose.pt"):
        self.video_path = video_path
        self.output_path = output_path
        self.mode = mode 
        self.num_workers = num_workers
        self.model_path = model_path
        
        # Буферы для producer-consumer
        self.input_buffer = queue.Queue(maxsize=10)
        self.output_buffer = queue.Queue()
        
        # Для восстановления порядка кадров
        self.frame_results: Dict[int, np.ndarray] = {}
        self.next_expected_id = 0
        
        # Сигналы управления
        self.stop_workers = threading.Event()
        self.workers: List[WorkerThread] = []
        self.fps = 25.0  
    
    def _read_video_info(self, cap: cv2.VideoCapture) -> tuple:
        width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
        height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
        fps = cap.get(cv2.CAP_PROP_FPS)
        self.fps = fps
        return width, height, fps
    
    def _producer(self, cap: cv2.VideoCapture):
        frame_id = 0
        while not self.stop_workers.is_set():
            ret, frame = cap.read()
            if not ret:
                break
            
            frame_data = FrameData(frame_id=frame_id, image=frame)
            self.input_buffer.put(frame_data)
            frame_id += 1
        
        # Сигнал завершения для воркеров
        for _ in range(self.num_workers):
            self.input_buffer.put(None)
    
    def _consumer(self, writer: cv2.VideoWriter):
        total_processed = 0
        
        while True:
            try:
                frame_data: FrameData = self.output_buffer.get(timeout=1.0)
            except queue.Empty:
                # Проверяем, не завершились ли все воркеры
                if all(not w.is_alive() for w in self.workers) and self.input_buffer.empty():
                    break
                continue
            
            # Сохраняем результат по ID
            self.frame_results[frame_data.frame_id] = frame_data.result
            
            # Записываем кадры в правильном порядке
            while self.next_expected_id in self.frame_results:
                result_frame = self.frame_results.pop(self.next_expected_id)
                writer.write(result_frame)
                self.next_expected_id += 1
                total_processed += 1
            
            self.output_buffer.task_done()

            if frame_data.frame_id == self.next_expected_id - 1 and not self.frame_results:
                if self.input_buffer.empty() and all(q.empty() for q in [self.input_buffer, self.output_buffer]):
                    break
        
        return total_processed
    
    def _single_thread_inference(self, cap: cv2.VideoCapture, writer: cv2.VideoWriter) -> float:
        model = YOLO(self.model_path)
        start_time = time.time()
        frame_id = 0
        
        while True:
            ret, frame = cap.read()
            if not ret:
                break
            
            results = model.predict(frame, verbose=False)
            if results and results[0].keypoints is not None:
                result_frame = results[0].plot()
            else:
                result_frame = frame.copy()
            
            writer.write(result_frame)
            frame_id += 1
        
        return time.time() - start_time
    
    def _multi_thread_inference(self, cap: cv2.VideoCapture, writer: cv2.VideoWriter) -> float:
        start_time = time.time()
        
        # Запуск воркеров
        for _ in range(self.num_workers):
            worker = WorkerThread(
                input_queue=self.input_buffer,
                output_queue=self.output_buffer,
                model_path=self.model_path,
                stop_event=self.stop_workers
            )
            worker.start()
            self.workers.append(worker)
        
        # Producer и Consumer в отдельных потоках
        producer_thread = threading.Thread(target=self._producer, args=(cap,), daemon=True)
        consumer_thread = threading.Thread(target=self._consumer, args=(writer,), daemon=True)
        
        producer_thread.start()
        consumer_thread.start()
        
        # Ожидание завершения
        producer_thread.join()
        for worker in self.workers:
            worker.join()
        consumer_thread.join()
        
        return time.time() - start_time
    
    def run(self) -> Dict[str, float]:
        results = {}
        
        with ResourceRAII("capture", source=self.video_path) as cap:
            width, height, fps = self._read_video_info(cap)
            
            with ResourceRAII("writer", filename=self.output_path, fps=fps, width=width, height=height) as writer:
                if self.mode == "single":
                    exec_time = self._single_thread_inference(cap, writer)
                else:
                    exec_time = self._multi_thread_inference(cap, writer)
                
                results["execution_time"] = exec_time
                results["fps"] = fps
        return results


def find_optimal_threads(video_path: str, model_path: str = "yolov8s-pose.pt", 
                         max_workers: int = 8, num_iterations: int = 10) -> int:
    print(f"\nПодбор оптимального числа потоков ({num_iterations} итераций)")
    best_avg_time = float('inf')
    best_workers = 1
    all_results = []

    for n_workers in range(1, min(max_workers + 1, 9)):
        iteration_times = []
        print(f"  Тестируем {n_workers} поток(ов)...")
        
        for i in range(num_iterations):
            temp_output = f"_bench_temp_{n_workers}_{i}.mp4"
            pipeline = PoseInferencePipeline(
                video_path=video_path,
                output_path=temp_output,
                mode="multi",
                num_workers=n_workers,
                model_path=model_path
            )
            
            start = time.time()
            pipeline.run()
            elapsed = time.time() - start
            iteration_times.append(elapsed)
            
            if os.path.exists(temp_output):
                try:
                    os.remove(temp_output)
                except PermissionError:
                    gc.collect()
                    time.sleep(0.2)
                    if os.path.exists(temp_output):
                        os.remove(temp_output)
            
            print(f"    [{i+1}/{num_iterations}] {elapsed:.2f}c", end='\r')
        print() 

        # Статистика
        avg_time = statistics.mean(iteration_times)
        std_dev = statistics.stdev(iteration_times) if len(iteration_times) > 1 else 0.0
        
        all_results.append((n_workers, avg_time, std_dev))
        
        if avg_time < best_avg_time:
            best_avg_time = avg_time
            best_workers = n_workers

    print("\nСводка:")
    print(f"  {'Потоков':<8} | {'Среднее время':<15} | {'Разброс (±)':<10} | Статус")
    for w, avg, std in all_results:
        marker = "ЛУЧШИЙ" if w == best_workers else "  "
        print(f"  {w:<8} | {avg:<15.2f} | {std:<10.2f} | {marker}")

    print(f"\nОптимальное число потоков: {best_workers} (среднее время: {best_avg_time:.2f}c)")
    return best_workers

def realtime_camera_demo(model_path: str = "yolov8s-pose.pt", num_workers: int = 4):
    print("Нажмите 'q' для выхода")
    
    cap = cv2.VideoCapture(0)
    cap.set(cv2.CAP_PROP_FRAME_WIDTH, 640)
    cap.set(cv2.CAP_PROP_FRAME_HEIGHT, 480)
    
    if not cap.isOpened():
        raise RuntimeError("Не удалось открыть камеру")
    
    input_q = queue.Queue(maxsize=5)
    output_q = queue.Queue(maxsize=10)
    
    stop_event = threading.Event()
    capture_done = threading.Event()
    
    workers = []
    for i in range(num_workers):
        w = WorkerThread(input_q, output_q, model_path, stop_event)
        w.start()
        workers.append(w)
    
    def capture_thread():
        frame_id = 0
        try:
            while not stop_event.is_set():
                ret, frame = cap.read()
                if not ret:
                    break
                try:
                    input_q.put(FrameData(frame_id=frame_id, image=frame), timeout=0.1)
                    frame_id += 1
                except queue.Full:
                    continue
        finally:
            capture_done.set()
    
    capture_t = threading.Thread(target=capture_thread, daemon=True)
    capture_t.start()
    
    display_buffer: Dict[int, np.ndarray] = {}
    next_display_id = 0
    start_time = time.time()
    frame_count = 0
    
    try:
        while not stop_event.is_set():
            try:
                fd: FrameData = output_q.get(timeout=0.05)  
            except queue.Empty:
                pass 
            else:
                display_buffer[fd.frame_id] = fd.result
                output_q.task_done()
            
            while next_display_id in display_buffer:
                frame = display_buffer.pop(next_display_id)
                cv2.imshow('123', frame)
                next_display_id += 1
                frame_count += 1
            
            key = cv2.waitKey(1) & 0xFF
            if key == ord('q') or key == 27: 
                print("\nПолучен сигнал выхода...")
                stop_event.set()
                break
            
            if frame_count % 30 == 0 and frame_count > 0:
                elapsed = time.time() - start_time
                current_fps = frame_count / elapsed
                print(f"\rReal-time FPS: {current_fps:.1f}", end='', flush=True)
    
    except KeyboardInterrupt:
        print("\nПрервано пользователем")
        stop_event.set()
    
    finally:
        capture_done.wait(timeout=1.0)
        
        stop_event.set()
        
        while not input_q.empty():
            try:
                input_q.get_nowait()
                input_q.task_done()
            except queue.Empty:
                break
        for _ in workers:
            input_q.put(None) 
        
        for i, w in enumerate(workers):
            w.join(timeout=2.0) 
            if w.is_alive():
                print(f"Воркер {i} не завершился, пропускаем...")
        
        while not output_q.empty():
            try:
                output_q.get_nowait()
                output_q.task_done()
            except queue.Empty:
                break
        cap.release()
        cv2.destroyAllWindows()


def main():
    parser = argparse.ArgumentParser(
        description="YOLOv8s-pose inference с поддержкой многопоточности",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Примеры использования:
  # Однопоточный режим
  python main.py --video video.mp4 --mode single --output result.mp4
  
  # Многопоточный режим с 4 потоками
  python main.py --video video.mp4 --mode multi --workers 4 --output result.mp4
  
  # Подбор оптимального числа потоков
  python main.py --video video.mp4 --mode auto --output result.mp4
  
  # Real-time с камеры
  python main.py --camera --workers 4
        """
    )
    
    parser.add_argument('--video', type=str)
    parser.add_argument('--mode', type=str, default='single', 
                       choices=['single', 'multi', 'auto'])
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--output', type=str, required=True,
                       help='Имя выходного видеофайла')
    parser.add_argument('--model', type=str, default='yolov8s-pose.pt')
    parser.add_argument('--camera', action='store_true')
    
    args = parser.parse_args()
    
    if args.camera:
        realtime_camera_demo(args.model, args.workers)
        return
    
    if not args.video:
        parser.error("Требуется указать --video или использовать --camera")
    
    if args.mode == 'auto':
        optimal = find_optimal_threads(args.video, args.model)
        args.workers = optimal
        args.mode = 'multi'
    
    print(f" Запуск: {args.mode} режим, {args.workers} потоков")
    pipeline = PoseInferencePipeline(
        video_path=args.video,
        output_path=args.output,
        mode=args.mode,
        num_workers=args.workers,
        model_path=args.model
    )
    
    start_total = time.time()
    stats = pipeline.run()
    total_time = time.time() - start_total

    print(f"\nРезультаты:")
    print(f"   Время обработки: {stats['execution_time']:.2f} секунд")
    print(f"   FPS видео: {stats['fps']:.1f}")
    print(f"   Общее время: {total_time:.2f} секунд")
    
    if args.mode == 'multi':
        print(f"   Потоков: {args.workers}")

if __name__ == "__main__":
    main()