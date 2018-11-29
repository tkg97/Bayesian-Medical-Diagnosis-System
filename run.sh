#!/bin/bash
g++ diagnoser.cpp -std=c++1y -o myoutput.out
./myoutput.out alarm1.bif records1.dat
g++ Format_Checker.cpp -std=c++1y -o checker.out
./checker.out
rm solved_alarm.bif myoutput.out checker.out