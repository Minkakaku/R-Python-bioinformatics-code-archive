cp /mnt/NexData_file/home/zwl/Project/personal/demo/SCENIC/new_code/report_scenic.rmd $1
cd $1
Rscript -e "rmarkdown::render('./report_scenic.rmd',output_file = 'report.html')"

