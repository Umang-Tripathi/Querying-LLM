// Do NOT add any other includes
#pragma once
#include <string> 
#include <vector>
#include <iostream>
#include "Node.h"

using namespace std;

class SearchEngine {
private:
    // You can add attributes/helper functions here
    vector<string> text_data;
    vector<Node* > sentence_info;
    vector<Node* > head_of_search;
public: 
    /* Please do not touch the attributes and 
    functions within the guard lines placed below  */
    /* ------------------------------------------- */
    SearchEngine();

    ~SearchEngine();

    void insert_sentence(int book_code, int page, int paragraph, int sentence_no, string sentence);

    Node* search(string pattern, int& n_matches);

    /* -----------------------------------------*/
};
