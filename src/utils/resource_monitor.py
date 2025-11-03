import time
import psutil
import subprocess

class ResourceMonitor:
    def run_with_monitoring(self, command):
        start_time = time.time()
        initial_memory = psutil.virtual_memory().used / (1024 ** 2)
        
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        peak_memory = initial_memory
        try:
            ps_process = psutil.Process(process.pid)
            while process.poll() is None:
                try:
                    mem_info = ps_process.memory_info()
                    current_mem = mem_info.rss / (1024 ** 2)
                    peak_memory = max(peak_memory, current_mem)
                except:
                    break
                time.sleep(0.1)
        except:
            peak_memory = 0
        
        stdout, stderr = process.communicate()
        end_time = time.time()
        
        return {
            'runtime': end_time - start_time,
            'peak_memory_mb': max(0, peak_memory - initial_memory),
            'return_code': process.returncode
        }
