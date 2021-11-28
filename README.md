Компиляция
==========

1. Перейти в папку проекта:

    `cd path/to/project`

2. Выполнить команду:

    `g++ main.cpp functions.cpp -fopenmp -o main.exe`

3. Фиксация изменений:
    `git add <file 1> <file 2> ...`
    `git commit -m "text of commit"`
    `git push origin master`
483 

254

g++ mpi_test.cpp -L"C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" -I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include" -o mpi_test.exe -lmsmpi
