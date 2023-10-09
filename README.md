# merge VCF
to merge multiple VCF

```shell
# usage:
python merge_VCF.py -h
```

```shell
# e.g., 
# same sites, different samples
python merge_VCF.py --filelist merge.list --merge sample --out out.vcf.gz
# same samples, different sites    
python merge_VCF.py --filelist merge.list --merge sites --out out.vcf.gz    
```

Note:    
merge.list is a file including a list of /path/to/vcf.gz    
developed in python3

---
By: Yuwen Pan  
Contact: [panyuwen.x@gmail.com](mailto:panyuwen.x@gmail.com)    

