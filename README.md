# similarity-constrained-TRCA (scTRCA)
**A modified task-related component analysis (TRCA) method was proposed for SSVEP detection, and the proposed method showed higher performance than traditional TRCA and other TRCA-based methods.**  
&nbsp;  
More details please refer to our paper: *Q. Sun, et al., “Similarity-constrained task-related component analysis for enhancing SSVEP detection,” J Neural Eng, vol. 18, no. 4, Jun 4, 2021.*
## Brief Description
- Demo code of scTRCA in MATLAB;
- Benchmark dataset from Tsinghua University was used for validation [1]；   
- Filter bank technique was adopted [2];   
- TRCA TRCA-R, and msTRCA were used for comparison [3-5].
## References
[1] Y. Wang, et al., “A benchmark dataset for SSVEP-based brain-computer interfaces,” IEEE Trans. Biomed. Eng., vol. 25, no. 10, pp. 1746-1752, Oct, 2017.  
[2] X. Chen, et al., “Filter bank canonical correlation analysis for implementing a high-speed SSVEP-based brain-computer interface,” J. Neural Eng., vol. 12, no. 4, pp. 046008, Aug, 2015.  
[3] M. Nakanishi, et al., “Enhancing detection of ssveps for a high-speed brain speller using task-related component analysis,” IEEE Trans. Biomed. Eng., vol. 65, no. 1, pp. 104-112, Jan, 2018.  
[4] C. M. Wong, et al., “Spatial filtering in SSVEP-based BCIs: Unified framework and new improvements,” IEEE Trans. Biomed. Eng., vol. 67, no. 11, pp. 3057-3072, Feb 21, 2020.  
[5] C. M. Wong, et al., “Learning across multi-stimulus enhances target recognition methods in SSVEP-based BCIs,” J. Neural Eng., vol. 17, no. 1, pp. 016026, Jan 6, 2020.  
## Acknowledgement
Thank Dr. Chi Man Wong very much for providing the code of TRCA-R (https://github.com/edwin465/SSVEP-TRCA-R).    
This demo was modified based on his data processing framework.
