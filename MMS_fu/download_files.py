import sys
import os
import time
import requests
from concurrent.futures import ThreadPoolExecutor

def download_chunk(args):
    url, start_byte, end_byte, local_file_path, file_name, progress_callback = args
    headers = {'Range': f'bytes={start_byte}-{end_byte}'}

    with requests.get(url, headers=headers, stream=True) as response:
        with open(local_file_path, 'r+b') as file:
            file.seek(start_byte)
            for chunk in response.iter_content(chunk_size=1024*10):
                if chunk:
                    file.write(chunk)
                    progress_callback(len(chunk))

def write_files(file_url, local_file_path, file_name, num_threads):
    response = requests.head(file_url)
    content_size = int(response.headers['content-length'])
    chunk_size = content_size // num_threads

    with open(local_file_path, 'wb') as file:
        file.truncate(content_size)

    def progress_callback(chunk_size):
        nonlocal current_size
        current_size += chunk_size
        elapsed_time = time.time() - start_time
        print(f'\rDownloaded: [{round(current_size / 1024 / 1024, 2)} MB]'
              f'[{round(float(current_size / content_size) * 100, 2)}%]'
              f'Elapsed Time: {round(elapsed_time, 2)} seconds', end="")

    start_time = time.time()
    current_size = 0
    futures = []
        
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for i in range(num_threads):
            start_byte = i * chunk_size
            end_byte = start_byte + chunk_size - 1 if i < num_threads - 1 else content_size - 1
            args = (file_url, start_byte, end_byte, local_file_path, file_name, progress_callback)
            futures.append(executor.submit(download_chunk, args))
    executor.shutdown(wait=True)
    
    for future in futures:
        future.result()

if __name__ == "__main__":
    file_url = sys.argv[1]
    local_file_path = sys.argv[2]
    file_name = sys.argv[3]
    num_threads = int(sys.argv[4])

    write_files(file_url, local_file_path, file_name, num_threads)
