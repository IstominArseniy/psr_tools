from tqdm.auto import tqdm
import time
from alive_progress import alive_bar
from alive_progress.styles import showtime

with alive_bar(1000, bar='scuba') as bar:  # your expected total
    for i in range(1000):        # the original loop
        time.sleep(1e-2)
        bar()         