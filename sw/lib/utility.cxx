#include<iostream>

#include "utility.h"

using namespace std;

void printProgress(double progress){
    int bar = 70;
    cout<<"[";
    int pos = bar*progress;
    for(int i = 0; i < bar; i++){
      if(i < pos) cout <<"=";
      else if(i == pos)cout  <<">";
      else cout <<" ";
    }
    cout <<"]  " <<int(progress*100) <<"%\r";
    cout.flush();
}

void printEnd(){
  int bar = 70;
  cout<<"[";
  int pos = bar;
  for(int i = 0; i < bar; i++){
    if(i < pos) cout <<"=";
    else if(i == pos)cout  <<">";
    else cout <<" ";
  }
  cout <<"] " <<int(100) <<"\n";
}
