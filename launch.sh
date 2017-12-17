#program launcher
cp out out_prev
ConfigFile=./Config/operational.cfg
main.exe ${ConfigFile} 1>out 2>err
