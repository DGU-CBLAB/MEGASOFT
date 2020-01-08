b2 variant=release  link=static runtime-link=static threading=multi address-model=64 -j4 stage

g++ -pthread -I ..\..\boost_1_72_0\ MetaSnp.cpp MetaSnp.h MetaSoft.cpp ..\..\boost_1_72_0\stage\lib\libboost_system-vc142-mt-x64-1_72.lib ..\..\boost_1_72_0\stage\lib\libboost_program_options-vc142-mt-x64-1_72.lib ..\..\boost_1_72_0\stage\lib\libboost_thread-vc142-mt-x64-1_72.lib -w -o MEGASOFT.exe