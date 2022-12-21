//
//  Progress.h
//  Tile
//
//  Created by T. Suzuki on 2013/08/18.
//  Copyright (c) 2013年 T. Suzuki. All rights reserved.
//

#ifndef __Tile__Progress__
#define __Tile__Progress__

#include <iostream>

using namespace std;

#define DIR_K (0x0001)      // 001
#define DIR_J (0x0002)      // 010
#define DIR_I (0x0004)      // 100
#define DONE  (0x0007)      // 111
#define NYET  (0x0000)      // 000

// Triplet Class
// This is the element of the task queue
class Triplet {
private:
	int i_;
	int j_;
	int k_;
	
public:
	Triplet() { i_=0; j_=0; k_=0; };
	Triplet(const int i, const int j, const int k) { i_ = i; j_ = j; k_ = k; };
	void setI(const int i) { i_ = i; };
	void setJ(const int j) { j_ = j; };
	void setK(const int k) { k_ = k; };
	void setIJK(const int i, const int j, const int k) { i_ = i; j_ = j; k_ = k; };
	void show() { cout << "(" << getI() << "," << getJ() << "," << getK() << ")"; };
	int getI() { return i_; };
	int getJ() { return j_; };
	int getK() { return k_; };
};

// Progress Table Class
class Progress_Table {
	
private:
	unsigned char* top_;     // pointer to the progress table
	int m_;     // number of lows
	int n_;     // number of columns
	int k_;     // depth
	
public:
	Progress_Table( const int M=1, const int N=1, const int K=1);
	void setIJK( const int, const int, const int, const unsigned char );
	unsigned char getIJK( const int, const int, const int );
	void check_waitIJK( const int, const int, const int );
	void Init();
	void Show();
};

#endif /* defined(__Tile__Progress__) */
