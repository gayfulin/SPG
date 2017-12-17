#program launcher
cp test_out test_out_prev
ConfigFile=./Config/testrun.cfg

main.exe ${ConfigFile} 1>test_out 2>test_err
