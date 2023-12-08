// Do NOT add any other includes
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#pragma once
using namespace std;
struct TrieNode {
    char key;
    long long freq;
    bool is_end;
    vector<TrieNode *> children;
    TrieNode(char key, long long freq, bool is_end)
        : key{key}, freq{freq}, is_end{is_end} {}
    ~TrieNode() {
        for (int i = 0; i < children.size(); i++) {
            TrieNode *child = children[i];
            delete child;
        }
    }
};

class Trie {
   public:
    Trie() { root = nullptr; }
    ~Trie() { delete root; }
    TrieNode *root;
    void increase_freq(string word) {
        if (!root) {
            // FIXME: Using - as zero length char, might cause issues
            // if "well formed sentences" have a hyphen in them.
            root = new TrieNode('-', 1, false);
        }
        increase_freq_help(root, word, 0);
    }
    void increase_freq_help(TrieNode *cur, string &word, int index) {
        TrieNode *child = nullptr;
        for (int i = 0; i < cur->children.size(); i++) {
            auto elem = cur->children[i];
            if (elem->key == word[index]) {
                child = elem;
                break;
            }
        }
        if (!child) {
            child = new TrieNode(word[index], 0, false);
            if (index == word.length() - 1) {
                child->is_end = true;
                child->freq = 1;
            } else {
                increase_freq_help(child, word, ++index);
            }
            cur->children.push_back(child);
        } else {
            if (index == word.length() - 1) {
                child->is_end = true;
                child->freq++;
            } else {
                increase_freq_help(child, word, ++index);
            }
        }
    }
    void set_freq(string word, long long freq) {
        if (!root) {
            // FIXME: Using - as zero length char, might cause issues
            // if "well formed sentences" have a hyphen in them.
            root = new TrieNode('-', 1, false);
        }
        set_freq_help(root, word, 0, freq);
    }
    void set_freq_help(TrieNode *cur, string word, int index, long long freq) {
        TrieNode *child = nullptr;
        for (int i = 0; i < cur->children.size(); i++) {
            auto elem = cur->children[i];
            if (elem->key == word[index]) {
                child = elem;
                break;
            }
        }
        if (!child) {
            child = new TrieNode(word[index], 0, false);
            if (index == word.length() - 1) {
                child->is_end = true;
                child->freq = freq;
            } else {
                set_freq_help(child, word, ++index, freq);
            }
            cur->children.push_back(child);
        } else {
            if (index == word.length() - 1) {
                child->is_end = true;
                child->freq = freq;
            } else {
                set_freq_help(child, word, ++index, freq);
            }
        }
    }

    long long get_freq(string word) {
        if (!root) {
            return 0;
        }
        return get_freq_help(root, word, 0);
    }
    long long get_freq_help(TrieNode *cur, string &word, int index) {
        TrieNode *child = nullptr;
        for (int i = 0; i < cur->children.size(); i++) {
            TrieNode *elem = cur->children[i];
            if (elem->key == word[index]) {
                child = elem;
                break;
            }
        }
        if (!child) {
            return 0;
        } else {
            if (index == word.length() - 1) {
                if (child->is_end) {
                    return child->freq;
                } else {
                    return 0;
                }
            } else {
                return get_freq_help(child, word, ++index);
            }
        }
    }
};
class Dict {
private:
  // You can add attributes/helper functions here
  Trie freq_table;

public:
  /* Please do not touch the attributes and
  functions within the guard lines placed below  */
  /* ------------------------------------------- */
  Dict();
  vector<int> dataIn;

  ~Dict();

  void insert_sentence(int book_code, int page, int paragraph, int sentence_no,
                       string sentence);

  int get_word_count(string word);

  void dump_dictionary(string filename);

  /* -----------------------------------------*/
};

string to_lower(string word);
bool is_separator(char c);
