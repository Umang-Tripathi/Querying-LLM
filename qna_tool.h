#pragma once
#include "Node.h"
#include "dict.h"
#include <fstream>
#include <iostream>
using namespace std;

class QNA_tool {

private:
  // You are free to change the implementation of this function
  void query_llm(string filename, Node *root, int k, string API_KEY);
  // You can add attributes/helper functions here
  vector<vector<vector<Trie *>>> ParaTrie;

  vector<vector<vector<Trie *>>> StemmedParaTrie;                                  //stemmed trie

  Trie *StandardText = new Trie();

  Trie *StemmedStandardText = new Trie();                                          //stemmed trie
  
  Trie *corpus = new Trie();

  Trie *StemmedCorpus = new Trie();                                                //stemmed trie
public:
  /* Please do not touch the attributes and
  functions within the guard lines placed below  */
  /* ------------------------------------------- */

  QNA_tool();  // Constructor
  ~QNA_tool(); // Destructor

  void insert_sentence(int book_code, int page, int paragraph, int sentence_no,
                       string sentence);
  // This function is similar to the functions in dict and search
  // The corpus will be provided to you via this function
  // It will be called for each sentence in the corpus

  Node *get_top_k_para(string question, int k);
  // This function takes in a question, preprocess it
  // And returns a list of paragraphs which contain the question
  // In each Node, you must set: book_code, page, paragraph (other parameters
  // won't be checked)

  std::string get_paragraph(int book_code, int page, int paragraph);
  // Given the book_code, page number, and the paragraph number, returns the
  // string paragraph. Searches through the corpus.
  void query(string question, string filename);
  // This function takes in a question and a filename.
  // It should write the final answer to the specified filename.

  void query_llm(string filename, Node *root, int k, string API_KEY,
                 string question);
  /* -----------------------------------------*/
  /* Please do not touch the code above this line */

  // You can add attributes/helper functions here
};
