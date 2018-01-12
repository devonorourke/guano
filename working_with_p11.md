original file structure:
```
NHCS-oahuRat261-xx-VE-USA-2017-076-JF_S294_L001_R1_001.fastq.gz
```

use rename to remove front end:
```
rename NHCS-oahuRat oahuRat *.gz
```

then rename again to remove middle portion:
```
rename -- -xx-VE-USA-2017-076-JF_S*_L001_ - *.gz
```
