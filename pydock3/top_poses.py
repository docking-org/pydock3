import sys, os
import subprocess as sp
import gzip
import time
import multiprocessing

from pydock3.util import Script


def get_to_search(dock_results_dir_path):

    if os.path.isfile(dock_results_dir_path):
        with open(dock_results_dir_path, 'r') as f:
            for line in f:
                yield line.strip()
    elif os.path.isdir(dock_results_dir_path):
        with sp.Popen(["find", "-L", dock_results_dir_path, "-type", "f", "-name", "test.mol2.gz.*"], stdout=sp.PIPE) as f:
            for line in f.stdout:
                yield line.decode('utf-8').strip()
    else:
        print(f"Supplied dock results path {dock_results_dir_path} cannot be found!")
        sys.exit(1)


class MinHeap(object):

    def __init__(self, max_size=10000, comparator=lambda x, y:x < y):
        self.max_size = max_size
        self.comparator = comparator
        self.size = 0
        self.heap = [None for i in range(max_size+1)]

    def __swap(self, i, j):
        tmp = self.heap[i]
        self.heap[i] = self.heap[j]
        self.heap[j] = tmp

    # best case: O(1)
    # worst case: O(log2(n))
    def insert(self, element):
        self.size += 1
        self.heap[self.size] = element
        if self.size == 1:
            self.heap[0] = element
            return

        current = self.size
        parent = self.size // 2

        while self.comparator(self.heap[current], self.heap[parent]):
            self.__swap(current, parent)
            current = parent
            parent = current // 2

        self.heap[0] = self.heap[1]

    # best case: O(1)
    # worst case: O(log2(n))
    def remove_insert(self, element):
        popped = self.heap[1]
        self.heap[1] = element

        position = 1
        while position < (self.size // 2):
            left_child = 2 * position
            right_child = 2 * position + 1

            if self.comparator(self.heap[position], self.heap[right_child]) and self.comparator(self.heap[position], self.heap[left_child]):
                break

            if self.comparator(self.heap[left_child], self.heap[right_child]):
                self.__swap(position, left_child)
                position = left_child
            else:
                self.__swap(position, right_child)
                position = right_child

        self.heap[0] = self.heap[1]
        return popped

    def minvalue(self):
        return self.heap[1]


class BufferedLineReader(object):

    def __init__(self, file_name, buffer_size=1024*1024):
        self.buffer_size = buffer_size
        if '.gz' in file_name:
            try:
                self.file = gzip.open(file_name, 'rt')
            except:
                self.file = open(file_name, 'r')
        else:
            self.file = open(file_name, 'r')
        self.partial_line = None
        self.buffer = None
        self.__read_buffer()

    def __read_buffer(self):
        del self.buffer
        self.buffer = self.file.read(self.buffer_size)
        self.buffer_position = 0
        self.buffer_length = len(self.buffer)

    def read_line(self):
        start = self.buffer_position
        while self.buffer_position < self.buffer_length and self.buffer[self.buffer_position] != '\n':
            self.buffer_position += 1
        if self.buffer_position == self.buffer_length: # we reached the end of the buffer without finding a newline
            if not self.buffer:
                return ''  # avoid an infinite recursive loop if buffer is empty
            partial_line = self.buffer[start:self.buffer_position]
            self.__read_buffer()
            full_line = partial_line + self.read_line() # bet you didn't expect to see recursion in a read_line function!
            # if this function is recursing more than once then you're dealing with some monster lines (or a very small buffer)
            return full_line
        self.buffer_position += 1
        return self.buffer[start:self.buffer_position]


class Mol2Data(object):

    def __init__(self, data, total_energy, name):
        self.buffer = data.encode('utf-8')
        self.total_energy = total_energy
        self.name = name


def pose_data_producer(processing_queue, paths_enum):
    for poses_file_path in paths_enum:
        try:
            reader = BufferedLineReader(poses_file_path)
            buff = ""
            energy = 99999
            name = ""
            line = reader.read_line()
        except:
            print(f"Encountered error while trying to open {poses_file_path}. Skipping!")
            return

        error = 0
        processing_queue.put(0)

        while line:
            if line.startswith("##########                 Name"):
                if buff:
                    processing_queue.put(Mol2Data(buff, energy, name))
                    buff = ""
                    energy = 99999
                name = line.split(':')[1].strip()
            elif line.startswith("##########         Total Energy"):
                energy = float(line.split(':')[1].strip())

            buff += line
            try:
                line = reader.read_line()
            except:
                print(f"Encountered error while reading {poses_file_path}. Stopping read.")
                error = 1

        if error == 0:
            processing_queue.put(Mol2Data(buff, energy, name))

        processing_queue.put(poses_file_path)
    processing_queue.put(1)
    
    
def energy_is_greater_than_other_energy(m_1, m_2):
    return m_1.total_energy > m_2.total_energy
    

class TopPoses(Script):
    def __init__(self, dock_results_dir_path):
        #
        super().__init__()

        #
        self.dock_results_dir_path = dock_results_dir_path

    def run(self, output_file_path="top_poses.mol2.gz", max_heap_size=10000):
        heap = MinHeap(max_heap_size, comparator=energy_is_greater_than_other_energy)
        processing_queue = multiprocessing.Queue(maxsize=50)
        pose_data_producer_process = multiprocessing.Process(target=pose_data_producer, args=(processing_queue, get_to_search(self.dock_results_dir_path)))

        pose_data_producer_process.start()

        i = 0
        hit_max = False
        while True:
            try:
                pose = processing_queue.get(True, 10)
            except:
                break

            # int signals start the clock
            if type(pose) == int:
                if pose == 0:  # timer start
                    i = 0
                    start_time_u = time.clock_gettime(time.CLOCK_THREAD_CPUTIME_ID)
                    start_time_r = time.time()
                if pose == 1: # end processing
                    break
                continue

            # str signals to stop the clock, and also print the file_name (which will be the item)
            elif type(pose) == str:
                end_time_u = time.clock_gettime(time.CLOCK_THREAD_CPUTIME_ID)
                end_time_r = time.time()
                elapsed_u = end_time_u - start_time_u
                elapsed_r = end_time_r - start_time_r
                print(i, "poses read from", pose)
                print("time (user): {:15f}, time (real): {:15f}".format(elapsed_u, elapsed_r))
                print("pps  (user): {:15f}, pps  (real): {:15f}".format(i / elapsed_u, i / elapsed_r))
                continue

            if heap.size == max_heap_size:
                if not hit_max:
                    print("hit max poses for heap!")
                    hit_max = True

                if energy_is_greater_than_other_energy(pose, heap.minvalue()):
                    del pose
                    continue

                else:
                    # pops the minimum value and inserts the new value in the same operation
                    popped = heap.remove_insert(pose)
                    del popped
            else:
                heap.insert(pose)

            i += 1

        print("done processing!")

        data_index = sorted([(i, heap.heap[i].total_energy) for i in range(heap.size)], key=lambda m: m[1])
        with gzip.open(output_file_path, 'w') as f:
            for idx, out_energy in data_index:
                # gzip data can simply be concatenated together, so this works fine.
                f.write(heap.heap[idx].buffer)
