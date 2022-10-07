wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz;
tar -zxvf sratoolkit.3.0.0-ubuntu64.tar.gz;
cd sratoolkit.3.0.0-ubuntu64/bin;
vdb-config -i;
sudo apt install sra-toolkit;
