// Do NOT add any other includes
#include "dict.h"
Dict::Dict() {
  // Implement your function here
}

Dict::~Dict() {
  // freq_table's destructor is called automatically
}

string to_lower(string word) {
  string lowerWord = "";
  for (int i = 0; i < word.length(); i++) {
    lowerWord += tolower(word[i]);
  }
  return lowerWord;
}
bool is_separator(char c) {
  switch (c) {
  case '.':
  case ',':
  case '-':
  case ':':
  case '\"':
  case '\'':
  case '!':
  case '(':
  case ')':
  case '?':
  case '[':
  case ']':
  case ';':
  case '@':
  case ' ':
    return true;
  default:
    return false;
  }
}
void Dict::insert_sentence(int book_code, int page, int paragraph,
                           int sentence_no, string sentence) {
  // FIXME: Using @Umang-Tripathi's old parser
  // Need to add support for punctuation and Viraj's "well-formed sentences"
  string word = "";
  for (int i = 0; i < sentence.size(); i++) {
    if (is_separator(sentence[i])) {
      if (word != "") {
        // cout << to_lower(word) << endl;
        this->freq_table.increase_freq(to_lower(word));
      }
      word = "";
    } else {
      word += sentence[i];
    }
  }
  if (word != "") {
    this->freq_table.increase_freq(to_lower(word));
  }
  return;
}
struct freq_record {
  string word;
  int freq;
  string get_printable() { return word + ", " + to_string(freq) + "\n"; }
};

int Dict::get_word_count(string word) {
  return this->freq_table.get_freq(to_lower(word));
}
void dfs_find_words(TrieNode *cur, string &prefix, vector<freq_record> &finds) {
  if (!cur) {
    return;
  }
  string new_prefix;
  // FIXME: Choose different character for freq
  if (cur->key != '-') {
    new_prefix = prefix + cur->key;
  } else {
    new_prefix = prefix;
  }
  if (cur->is_end) {
    finds.push_back(freq_record{new_prefix, (int)cur->freq});
  }
  for (int i = 0; i < cur->children.size(); i++) {
    auto child = cur->children[i];
    dfs_find_words(child, new_prefix, finds);
  }
}

void Dict::dump_dictionary(string filename) {
  // Implement your function here
  // Current idea is simple DFS over the entire trie.
  // Discuss with @Umang-Tripathi
  vector<freq_record> finds;
  string prefix = "";
  dfs_find_words(this->freq_table.root, prefix, finds);
  ofstream write_file;
  write_file.open(filename);
  // cout << filename << endl;
  for (int i = 0; i < finds.size(); i++) {
    freq_record entry = finds[i];
    write_file << entry.get_printable();
  }
  write_file.close();
  return;
}



