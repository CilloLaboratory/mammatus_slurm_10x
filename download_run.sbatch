#!/bin/bash
## NOTE: this could be simplified and updated using xargs based on https://www.baeldung.com/linux/wget-parallel-downloading
## NOTE: expects the directory to hold the fastq files to be provided as an environmental variable
## FASTQ_DIR, name of directory to hold FASTQ files
DOWNLINKS=($(ls | grep -Z "Download_Links"))

cut -d "," -f3 $DOWNLINKS | sed 1d > expected_md5.txt

cut -d "," -f2 $DOWNLINKS | sed 1d > download_links.txt

cut -d "," -f1 $DOWNLINKS | sed 1d > file_names.txt

FILE_NAME=($(while IFS= read -r line
  do
    echo "$line"
  done < file_names.txt
  ))

FILE=($(while IFS= read -r line
  do
    echo "$line"
  done < download_links.txt
  ))

rm file_names.txt download_links.txt

LENGTH=${#FILE[@]}

mkdir -p ${FASTQ_DIR}
cd ${FASTQ_DIR}

for (( i=0; i<${LENGTH}; i++ ));
  do
    echo "Downloading
    ${FILE[i]}
    "
    wget ${FILE[i]}
    if [ $i == 0 ]
      then
        md5sum ${FILE_NAME[i]} | cut -d" " -f1 > downloaded_md5.txt
      else
        md5sum ${FILE_NAME[i]} | cut -d" " -f1 >> downloaded_md5.txt
    fi
  done
