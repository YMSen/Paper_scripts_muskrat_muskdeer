//---------------------------------------------------------------------------
#include "gcl/GString.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "gcl/GBase.h"
#include <stdarg.h>
//---------------------------------------------------------------------------

GString::Data GString::null_data;

//=========================================

GString::Data * GString::new_data(int length) {
//static method to return a new Data object (allocate length)
//content is undefined, but it's null terminated
    if (length > 0) {
        Data* data;
        GMALLOC(data, sizeof(Data)+length);
        data->ref_count = 0;
        data->length = length;
        data->chars[length] = '\0';
        return data;
        }
    else
        return &null_data;
 }

GString::Data* GString::new_data(const char* str) {
//static method to return a new Data object (allocate length)
//as a copy of a given string
 if (str==NULL) return &null_data;
 int length=strlen(str);
 if (length > 0) {
        Data* data;
        GMALLOC(data, sizeof(Data)+length);
        strcpy(data->chars, str);
        data->ref_count = 0;
        data->length = length;
        data->chars[length] = '\0';
        return data;
        }
    else
        return &null_data;
 }
 
void GString::replace_data(int len) {

    if (len == my_data->length && my_data->ref_count <= 1)
        return;

    if (my_data != &null_data && --my_data->ref_count == 0)
        GFREE(my_data);

    if (len > 0) {
        //my_data = (Data *) malloc(sizeof(Data) + len);
        GMALLOC(my_data, sizeof(Data) + len);
        my_data->ref_count = 1;
        my_data->length = len;
        my_data->chars[len] = '\0';
    }
    else
        my_data = &null_data;
}

void GString::replace_data(Data *data) {
    if (my_data != &null_data && --my_data->ref_count == 0)
        GFREE(my_data);
    if (data != &null_data)
        data->ref_count++;
    my_data = data;
}

void GString::make_unique() {//make sure is not a reference to other string
    if (my_data->ref_count > 1) {
        Data *data = new_data(length());
        ::memcpy(data->chars, chars(), length());
        my_data->ref_count--;
        my_data = data;
        my_data->ref_count++;
    }
}

bool operator==(const char *s1, const GString& s2){
  if (s1==NULL) return s2.is_empty();
  return (strcmp(s1, s2.chars()) == 0);
  }

bool operator<(const char *s1, const GString& s2) {
  if (s1==NULL) return !s2.is_empty();
  return (strcmp(s1, s2.chars()) < 0);
  }

bool operator<=(const char *s1, const GString& s2){
 if (s1==NULL) return true;
 return (strcmp(s1, s2.chars()) <= 0);
 }

bool operator>(const char *s1, const GString& s2) {
  if (s1==NULL) return false;
 return (strcmp(s1, s2.chars()) > 0);
 }


GString::GString():my_data(&null_data) {
 fTokenDelimiter=NULL;
 fLastTokenStart=0;
 readbuf=NULL;
 }

GString::GString(const GString& s): my_data(&null_data){
 fTokenDelimiter=NULL;
 fLastTokenStart=0;
 readbuf=NULL;
 replace_data(s.my_data);
 }

GString::GString(const char *s): my_data(&null_data) {
  fTokenDelimiter=NULL;
  fLastTokenStart=0;
  readbuf=NULL;
  my_data=new_data(s);
  my_data->ref_count = 1;
 }

GString::GString(const int i): my_data(&null_data) {
 fTokenDelimiter=NULL;
 fLastTokenStart=0;
 readbuf=NULL;
 char buf[20];
 sprintf(buf,"%d",i);
 const int len = ::strlen(buf);
 replace_data(len);
 ::memcpy(chrs(), buf, len);
 }

GString::GString(const double f): my_data(&null_data) {
 fTokenDelimiter=NULL;
 fLastTokenStart=0;
 readbuf=NULL;
 char buf[20];
 sprintf(buf,"%f",f);
 const int len = ::strlen(buf);
 replace_data(len);
 ::memcpy(chrs(), buf, len);
 }

GString::GString(char c, int n): my_data(&null_data) {
  fTokenDelimiter=NULL;
  fLastTokenStart=0;
  readbuf=NULL;
  replace_data(n); ::memset(chrs(), c, n);
  }

GString::~GString() {  
  if (my_data != &null_data && --my_data->ref_count == 0)
             GFREE(my_data);
  GFREE(fTokenDelimiter);
  GFREE(readbuf);
  }

char& GString::operator[](int idx){
//returns reference to char (can be l-value)
  if (idx < 0) idx += length();
  if (idx < 0 || idx >= length()) invalid_index_error("operator[]");
  make_unique();  //because the user will probably modify this char!
  return chrs()[idx]; 
  }

char GString::operator[](int idx) const {
//returns char copy (cannot be l-value!)
  if (idx < 0) idx += length();
  if (idx < 0 || idx >= length()) invalid_index_error("operator[]");
  return chars()[idx];
  }

GString& GString::operator=(const GString& s) {
  make_unique(); //edit operation ahead
  replace_data(s.my_data); 
  return *this;
  }

GString& GString::operator=(const char *s) {
  make_unique(); //edit operation ahead
  if (s==NULL) {
    replace_data(0);
    return *this;
    }
  const int len = ::strlen(s); replace_data(len);
  ::memcpy(chrs(), s, len); 
  return *this;
  }

GString& GString::operator=(const double f) {
 make_unique(); //edit operation ahead
 char buf[20];
 sprintf(buf,"%f",f);
 const int len = ::strlen(buf);
 replace_data(len);
 ::memcpy(chrs(), buf, len);
 return *this;
}

GString& GString::operator=(const int i) {
 make_unique(); //edit operation ahead
 char buf[20];
 sprintf(buf,"%d",i);
 const int len = ::strlen(buf);
 replace_data(len);
 ::memcpy(chrs(), buf, len);
 return *this;
}

bool GString::operator==(const GString& s) const {
  if (s.is_empty()) return is_empty();
  return (length() == s.length()) &&
    (memcmp(chars(), s.chars(), length()) == 0);
  }

bool GString::operator==(const char *s) const {
 if (s==NULL) return is_empty();
 return (strcmp(chars(), s) == 0);
 }

bool GString::operator<(const GString& s) const {
 if (s.is_empty()) return false;
 return (strcmp(chars(), s.chars()) < 0);
 }

bool GString::operator<(const char *s) const {
 if (s==NULL) return false;
 return (strcmp(chars(), s) < 0);
 }

bool GString::operator<=(const GString& s) const {
 if (s.is_empty()) return is_empty();
 return (strcmp(chars(), s.chars()) <= 0);
 }

bool GString::operator<=(const char *s) const {
 if (s==NULL) return is_empty();
 return (strcmp(chars(), s) <= 0);
 }

bool GString::operator>(const GString& s) const {
 if (s.is_empty()) return !is_empty();
 return (strcmp(chars(), s.chars()) > 0);
 }

bool GString::operator>(const char *s) const {
 if (s==NULL) return !is_empty();
 return (strcmp(chars(), s) > 0);
 }

bool GString::operator>=(const GString& s) const {
 if (s.is_empty()) return true;
 return (strcmp(chars(), s.chars()) >= 0);
 }

bool GString::operator>=(const char *s) const {
 if (s==NULL) return true;
 return (strcmp(chars(), s) >= 0);
 }

bool GString::operator!=(const GString& s) const {
  if (s.is_empty()) return !is_empty();
  return (length() != s.length()) ||
         (memcmp(chars(), s.chars(), length()) != 0);
  }

bool GString::operator!=(const char *s) const {
 if (s==NULL) return !is_empty();
 return (strcmp(chars(), s) != 0);
 }

GString& GString::operator+=(const GString& s) {
 return append((const char *)s);
 }

GString& GString::operator+=(const char* s) {
 return append(s);
 }

GString& GString::operator+=(const char c) {
 char buf[4];
 sprintf(buf,"%c",c);
 return append(buf);
 }

GString& GString::operator+=(const int i) {
 char buf[20];
 sprintf(buf,"%d",i);
 return append(buf);
 }


GString& GString::operator+=(const double f) {
 char buf[30];
 sprintf(buf,"%f",f);
 return append(buf);
 }
 
bool GString::is_empty() const {
  //return my_data == &null_data;
  return (length()==0);
  }

GString GString::copy() const {
 GString newstring(*this);
 return newstring;
 }

GString& GString::clear() {
  make_unique(); //edit operation ahead
  replace_data(0);
  return *this;
  }

int GString::index(const GString& s, int start_index) const {
 return index(s.chars(), start_index);
 }

bool GString::contains(const GString& s) const {
 return (index(s, 0) >= 0);
 }

bool GString::contains(const char *s) const {
 return (index(s, 0) >= 0);
 }

bool GString::contains(char c) const {
 return (index(c, 0) >= 0);
 }
GString& GString::format(const char *fmt,...) {
// Format as in sprintf
  make_unique(); //edit operation ahead
  char* buf;
  GMALLOC(buf, strlen(fmt)+1024);
  va_list arguments;
  va_start(arguments,fmt);
  //+1K buffer, should be enough for common expressions
  int len=vsprintf(buf,fmt,arguments);
  va_end(arguments);
  replace_data(len); //this also adds the '\0' at the end!
                     //and sets the right len
  ::memcpy(chrs(), buf, len);
  GFREE(buf);
  return *this;
  }

GString& GString::appendfmt(const char *fmt,...) {
// Format as in sprintf
  make_unique(); //edit operation ahead
  char* buf;
  GMALLOC(buf, strlen(fmt)+1024);
  va_list arguments;
  va_start(arguments,fmt);
  //+1K buffer, should be enough for common expressions
  vsprintf(buf,fmt,arguments);
  va_end(arguments);
  append(buf);
  GFREE(buf);
  return *this;
  }

GString& GString::trim(char c) {
 register int istart;
 register int iend;
 for (istart=0; istart<length() && chars()[istart]==c;istart++);
 if (istart==length()) {
       make_unique(); //edit operation ahead
       replace_data(0); //string was entirely trimmed
       return *this;
       }
 for (iend=length()-1; iend>istart && chars()[iend]==c;iend--);
 int newlen=iend-istart+1;
 if (newlen==length())  //nothing to trim
           return *this; 
 make_unique(); //edit operation ahead
 Data *data = new_data(newlen);
 ::memcpy(data->chars, &chars()[istart], newlen);
 replace_data(data);
 return *this;
 }

GString& GString::trim(char* c) {
 register int istart;
 register int iend;
 for (istart=0; istart<length() && strchr(c, chars()[istart])!=NULL ;istart++);
 if (istart==length()) {
        replace_data(0); //string was entirely trimmed
        return *this;
        }
 for (iend=length()-1; iend>istart && strchr(c, chars()[iend])!=NULL;iend--);
 int newlen=iend-istart+1;
 if (newlen==length())  //nothing to trim
           return *this; 
 make_unique(); //edit operation ahead
 Data *data = new_data(newlen);
 ::memcpy(data->chars, &chars()[istart], newlen);
 replace_data(data);
 return *this;
 }

GString& GString::trimR(char c) {
 //only trim the right end
 //register int istart;
 register int iend;
 for (iend=length()-1; iend>=0 && chars()[iend]==c;iend--);
 if (iend==-1) {
       replace_data(0); //string was entirely trimmed
       return *this;
       }
 int newlen=iend+1;
 if (newlen==length())  //nothing to trim
           return *this; 
 make_unique(); //edit operation ahead

 Data *data = new_data(newlen);
 ::memcpy(data->chars, chars(), newlen);
 replace_data(data);
 return *this;
 }

GString& GString::trimR(char* c) {
 register int iend;
 for (iend=length()-1; iend>=0 && strchr(c,chars()[iend])!=NULL;iend--);
 if (iend==-1) {
       replace_data(0); //string was entirely trimmed
       return *this;
       }
 int newlen=iend+1;
 if (newlen==length())  //nothing to trim
           return *this; 
 make_unique(); //edit operation ahead
 Data *data = new_data(newlen);
 ::memcpy(data->chars, chars(), newlen);
 replace_data(data);
 return *this;
 }

GString& GString::trimL(char c) {
 register int istart;
 for (istart=0; istart<length() && chars()[istart]==c;istart++);
 if (istart==length()) {
       replace_data(0); //string was entirely trimmed
       return *this;
       }
 int newlen=length()-istart;
 if (newlen==length())  //nothing to trim
           return *this; 
 make_unique(); //edit operation ahead
 Data *data = new_data(newlen);
 ::memcpy(data->chars, &chars()[istart], newlen);
 replace_data(data);
 return *this;
 }

GString& GString::trimL(char* c) {
 register int istart;
 for (istart=0; istart<length() && strchr(c,chars()[istart])!=NULL;istart++);
 if (istart==length()) {
       replace_data(0); //string was entirely trimmed
       return *this;
       }
 int newlen=length()-istart;
 if (newlen==length())  //nothing to trim
           return *this; 
 make_unique(); //edit operation ahead

 Data *data = new_data(newlen);
 ::memcpy(data->chars, &chars()[istart], newlen);
 replace_data(data);
 return *this;
 }

GString& GString::padR(int len, char c) {
 //actually means align right in len
 if (length()>=len) return *this; //no room for padding
 make_unique(); //edit operation ahead
 Data *data = new_data(len);
 ::memset(data->chars,c,len-length());
 ::memcpy(&data->chars[len-length()], chars(), length());
 replace_data(data);
 return *this;
 }

GString& GString::padL(int len, char c) { //align left the string
 if (length()>=len) return *this; //no room for padding
 make_unique(); //edit operation ahead
 Data *data = new_data(len);
 ::memcpy(data->chars, chars(), length());
 ::memset(&data->chars[length()],c,len-length());
 replace_data(data);
 return *this;
 }

GString& GString::padC(int len, char c) {
 if (length()>=len) return *this; //no room for padding
 make_unique(); //edit operation ahead
 int istart=(len-length())/2;
 Data *data = new_data(len);
 if (istart>0)
      ::memset(data->chars, c, istart);
 ::memcpy(&data->chars[istart], chars(), length());
 int iend=istart+length();
 if (iend<len)
      ::memset(&data->chars[iend],c,len-iend);
 replace_data(data);
 return *this;
 }

GString operator+(const char *s1, const GString& s2) {
    const int s1_length = ::strlen(s1);

    if (s1_length == 0)
        return s2;
    else {
        GString newstring;
        newstring.replace_data(s1_length + s2.length());
        ::memcpy(newstring.chrs(), s1, s1_length);
        ::memcpy(&(newstring.chrs())[s1_length], s2.chars(), s2.length());
        return newstring;
        }
}

//=========================================

GString GString::operator+(const GString& s) const {
    if (length() == 0)
        return s;
    else if (s.length() == 0)
        return *this;
    else {
        GString newstring;
        newstring.replace_data(length() + s.length());
        ::memcpy(newstring.chrs(), chars(), length());
        ::memcpy(&(newstring.chrs())[length()], s.chars(), s.length());
        return newstring;
        }
}

//=========================================

GString GString::operator+(const char *s) const {

    const int s_length = ::strlen(s);

    if (s_length == 0)
        return *this;
    else {
        GString newstring;
        newstring.replace_data(length() + s_length);
        ::memcpy(newstring.chrs(), chars(), length());
        ::memcpy(&(newstring.chrs())[length()], s, s_length);
        return newstring;
        }
}

GString GString::operator+(const int i) const {
    char buf[20];
    sprintf(buf, "%d", i);
    const int s_length = ::strlen(buf);
    GString newstring;
    newstring.replace_data(length() + s_length);
    ::memcpy(newstring.chrs(), chars(), length());
    ::memcpy(&(newstring.chrs())[length()], buf, s_length);
    return newstring;
}

GString GString::operator+(const char c) const {
    char buf[4];
    sprintf(buf, "%c", c);
    const int s_length = ::strlen(buf);
    GString newstring;
    newstring.replace_data(length() + s_length);
    ::memcpy(newstring.chrs(), chars(), length());
    ::memcpy(&(newstring.chrs())[length()], buf, s_length);
    return newstring;
}

GString GString::operator+(const double f) const {
    char buf[30];
    sprintf(buf, "%f", f);
    const int s_length = ::strlen(buf);
    GString newstring;
    newstring.replace_data(length() + s_length);
    ::memcpy(newstring.chrs(), chars(), length());
    ::memcpy(&(newstring.chrs())[length()], buf, s_length);
    return newstring;
}


//=========================================

bool GString::is_space() const {

    if (my_data == &null_data)
        return false;

    for (register const char *p = chars(); *p; p++)
        if (!isspace(*p))
            return false;

    return true;
}

//=========================================

GString GString::substr(int idx, int len) const {
    // A negative idx specifies an idx from the right of the string.
    if (idx < 0)
        idx += length();

    // A length of -1 specifies the rest of the string.
    if (len == -1 || len>length()-idx)
        len = length() - idx;
    
    if (idx<0 || idx>=length() || len<0 )
        invalid_args_error("substr()");

    GString newstring;
    newstring.replace_data(len);
    ::memcpy(newstring.chrs(), &chars()[idx], len);
    return newstring;
}


//transform: any character from 'from' is replaced with a coresponding
//char from 'to'

GString&  GString::tr(char *rfrom, char* rto) {
 if (length() == 0 || rfrom==NULL || strlen(rfrom)==0)
        return *this;
 unsigned int l=strlen(rfrom);       
 if (rto!=NULL && strlen(rto)!=l)
      invalid_args_error("tr()");
 make_unique(); //edit operation ahead
 Data *data = new_data(length());
      
 if (rto==NULL) { //deletion case 
   char* s = my_data->chars;
   char* p;
   char* dest = data->chars;
   do {
      if ((p=strpbrk(s,rfrom))!=NULL) {
        memcpy(dest,s,p-s);
        dest+=p-s;
        s=p+1;
        }
       else { 
        strcpy(dest, s); 
        dest+=strlen(s);
        }
      } while (p!=NULL);
   (*dest)='\0';   
   }
  else { //char substitution case - easier!
   char* p;
   for (int i=0; i<length(); i++) {
    if ((p=strchr(rfrom, my_data->chars[i]))!=NULL) 
         my_data->chars[i]=rto[p-rfrom];
    }
   }
 data->length=strlen(data->chars);
 replace_data(data);
 return *this;
}


// search and replace all the occurences of a string with another string
// or just remove the given string (if replacement is NULL)
GString&  GString::replace(const char *rfrom, const char* rto) {
 if (length() == 0 || rfrom==NULL || strlen(rfrom)==0)
        return *this;
 unsigned int l=strlen(rfrom);
 unsigned int tl= (rto==NULL)?0:strlen(rto);
 make_unique(); //edit operation ahead
 char* p;
 char* dest;
 char* newdest=NULL;
 char* s = my_data->chars;
 if (tl!=l) { //reallocation
   if (tl>l) {  //possible enlargement
       GMALLOC(newdest, length()*(tl-l+1)+1);
       }
      else  {//delete or replace with a shorter string
       GMALLOC(newdest, length() + 1);
       }
      dest=newdest; 
      if (tl==0) {//deletion
           while ((p=strstr(s,rfrom))!=NULL) {
               //rfrom found at position p
                memcpy(dest,s,p-s);
                dest+=p-s;
                s+=p-s+l; //s positioned in string after rfrom
                }
           //no more occurences, copy the remaining string
           strcpy(dest, s);
          }
        else { //replace with another string
          while ((p=strstr(s,rfrom))!=NULL) {
              memcpy(dest,s,p-s); //copy up rto the match
              dest+=p-s;
              memcpy(dest,rto,tl); //put the replacement string
              dest+=tl;
              s+=p-s+l;
              }
          //not found any more, copy rto end of string
          strcpy(dest, s);
          }
       Data* data=new_data(newdest);
       replace_data(data);
       GFREE(newdest);
       }
  else { //inplace editing: no need rto reallocate
    while ((p=strstr(s,rfrom))!=NULL) {
        memcpy(p,rto,l);
        s+=p-s+l;
        }    
   }
 return *this;
}



GString&  GString::cut(int idx, int len) {

    if (len == 0)
        return *this;
    make_unique(); //edit operation ahead

    // A negative idx specifies an idx from the right of the string,
    // so the left part will be cut out
    if (idx < 0)
        idx += length();

    // A length of -1 specifies the rest of the string.
    if (len == -1)
        len = length() - idx;

    if (idx<0 || idx>=length() || len<0 || len>length()-idx)
        invalid_args_error("cut()");

    Data *data = new_data(length() - len);
    if (idx > 0)
        ::memcpy(data->chars, chars(), idx);
    ::strcpy(&data->chars[idx], &chars()[idx+len]);
    replace_data(data);

    return *this;
}

//=========================================

GString&  GString::paste(const GString& s, int idx, int len) {
    // A negative idx specifies an idx from the right of the string.
    if (idx < 0)
        idx += length();
    make_unique(); //edit operation ahead

    // A length of -1 specifies the rest of the string.
    if (len == -1)
        len = length() - idx;

    if (idx<0 || idx>=length() || len<0 || len>length()-idx)
        invalid_args_error("replace()");

    if (len == s.length() && my_data->ref_count == 1)
        ::memcpy(&chrs()[idx], s.chars(), len);
    else {
        Data *data = new_data(length() - len + s.length());
        if (idx > 0)
            ::memcpy(data->chars, chars(), idx);
        if (s.length() > 0)
            ::memcpy(&data->chars[idx], s.chars(), s.length());
        ::strcpy(&data->chars[idx+s.length()], &chars()[idx+len]);
        replace_data(data);
    }

    return *this;
}

//=========================================

GString& GString::paste(const char *s, int idx, int len) {

    // A negative idx specifies an idx from the right of the string.
    make_unique(); //edit operation ahead
    if (idx < 0)
        idx += length();

    // A length of -1 specifies the rest of the string.
    if (len == -1)
        len = length() - idx;

    if (idx<0 || idx>=length() || len<0 || len>length()-idx)
        invalid_args_error("replace()");

    const int s_length = ::strlen(s);

    if (len == s_length && my_data->ref_count == 1)
        ::memcpy(&chrs()[idx], s, len);
    else {
        Data *data = new_data(length() - len + s_length);
        if (idx > 0)
            ::memcpy(data->chars, chars(), idx);
        if (s_length > 0)
            ::memcpy(&data->chars[idx], s, s_length);
        ::strcpy(&data->chars[idx+s_length], &chars()[idx+len]);
        replace_data(data);
    }

    return *this;
}

//=========================================

GString& GString::insert(const GString& s, int idx) {
    make_unique(); //edit operation ahead

    // A negative idx specifies an idx from the right of the string.
    if (idx < 0)
        idx += length();

    if (idx < 0 || idx >= length())
        invalid_index_error("insert()");

    if (s.length() > 0) {
        Data *data = new_data(length() + s.length());
        if (idx > 0)
            ::memcpy(data->chars, chars(), idx);
        ::memcpy(&data->chars[idx], s.chars(), s.length());
        ::strcpy(&data->chars[idx+s.length()], &chars()[idx]);
        replace_data(data);
    }

    return *this;
}

//=========================================

GString& GString::insert(const char *s, int idx) {
    // A negative idx specifies an idx from the right of the string.
    make_unique(); //edit operation ahead
    if (idx < 0)
        idx += length();

    if (idx < 0 || idx >= length())
        invalid_index_error("insert()");

    const int s_length = ::strlen(s);

    if (s_length > 0) {
        Data *data = new_data(length() + s_length);
        if (idx > 0)
            ::memcpy(data->chars, chars(), idx);
        ::memcpy(&data->chars[idx], s, s_length);
        ::strcpy(&data->chars[idx+s_length], &chars()[idx]);
        replace_data(data);
    }

    return *this;
}
//=========================================

GString& GString::append(const char* s) {
  make_unique(); //edit operation ahead
  int len=::strlen(s);
  int newlength=len+my_data->length;
  if (newlength<=my_data->length) return *this;
  if (my_data->length==0) {
    replace_data(len);
    ::memcpy(my_data->chars, s, len);
    return *this;
   }
  //faster solution with realloc
  GREALLOC(my_data, sizeof(Data)+newlength);
  ::strcpy(&my_data->chars[my_data->length], s);
  my_data->length=newlength;
  my_data->chars[newlength]='\0';
  return *this;
}

GString& GString::append(const GString& s) {
 return append((const char *)s);
}


GString& GString::upper() {
  make_unique(); //edit operation ahead
  for (register char *p = chrs(); *p; p++)
            *p = (char) toupper(*p);

    return *this;
}

//=========================================

GString& GString::lower() {
    make_unique();

    for (register char *p = chrs(); *p; p++)
          *p = (char) tolower(*p);

    return *this;
}

//=========================================

int GString::index(const char *s, int start_index) const {
    // A negative index specifies an index from the right of the string.
    if (strlen(s)>(size_t)length()) return -1;
    if (start_index < 0)
        start_index += length();

    if (start_index < 0 || start_index >= length())
        invalid_index_error("index()");
    const char* idx = strstr(&chars()[start_index], s);
    if (!idx)
        return -1;
    else
        return idx - chars();
}

//=========================================

int GString::index(char c, int start_index) const {
    // A negative index specifies an index from the right of the string.
    if (length()==0) return -1;
    if (start_index < 0)
        start_index += length();
     
    if (start_index < 0 || start_index >= length())
        invalid_index_error("index()");


    if (c == '\0')
        return -1;
    const char *idx=(char *) ::memchr(&chars()[start_index], c,
                                         length()-start_index);
    if (!idx)
        return -1;
    else
        return idx - chars();
}

int GString::rindex(char c) const {   
    if (c == '\0' || length()==0)
        return -1;
    char* idx= rstrchr((char*)chars(), c);
    if (idx==NULL) return -1;
                else return idx-chars();
}

int GString::rindex(char* str) const {
    if (str==NULL || *str == '\0' || length()==0)
        return -1;
    char* idx= rstrfind((char*)chars(), str);
    if (idx==NULL) return -1;
                else return idx-chars();
}

GString GString::split(char* delim) {
           /* splits "this" in two parts, at the first (left) 
                 encounter of delim:
                 1st would stay in "this",
                 2nd part will be returned 
                 as a new string!
           */
 GString result;
 int i=index(delim);
 if (i>=0){
      result=substr(i+strlen(delim));
      cut(i);
      return result;
      }
 return result;
}

GString GString::split(char c) {
           /* splits "this" in two parts, at the first (left) 
                 encounter of delim:
                 1st would stay in "this",
                 2nd part will be returned 
                 as a new string!
           */
 GString result;
 int i=index(c);
 if (i>=0){
      result=substr(i+1);
      cut(i);
      return result;
      }
 return result;
}

GString GString::splitr(char* delim) {
 GString result;
 int i=rindex(delim);
 if (i>=0){
      result=substr(i+strlen(delim));
      cut(i);
      return result;
      }
 return result;
}

GString GString::splitr(char c) {
 GString result;
 int i=rindex(c);
 if (i>=0){
      result=substr(i+1);
      cut(i);
      return result;
      }
 return result;
}


void GString::startTokenize(const char* delimiter, enTokenizeMode tokenizemode) {
 GFREE(fTokenDelimiter);
 GMALLOC(fTokenDelimiter,strlen(delimiter)+1);
 strcpy(fTokenDelimiter, delimiter);
 fLastTokenStart=0;
 fTokenizeMode=tokenizemode;
}

bool GString::nextToken(GString& token) {
 if (fTokenDelimiter==NULL) {
    GError("GString:: no token delimiter; use StartTokenize first\n");
    }
 if (fLastTokenStart>=length()) {//no more
    GFREE(fTokenDelimiter);
    fLastTokenStart=0;
    return false;
    }
 int dlen=strlen(fTokenDelimiter);
 char* delpos=NULL; //delimiter position
 int tlen=0;
 if (fTokenizeMode==tkFullString) { //exact string as a delimiter
   delpos=(char*)strstr(chars()+fLastTokenStart,fTokenDelimiter);
   if (delpos==NULL) delpos=(char*)(chars()+length());
   //empty records may be returned
   if (chars()+fLastTokenStart == delpos) { //empty token
     fLastTokenStart=(delpos-chars())+dlen;
     token="";
     return true;
     }
    else {
     tlen=delpos-(chars()+fLastTokenStart);
     token.replace_data(tlen);
     ::memcpy(token.chrs(), &chars()[fLastTokenStart], tlen);
     fLastTokenStart=(delpos-chars())+dlen;
     return true;
     } 
   }
  else { //tkCharSet - any character is a delimiter
   //empty records are never returned !
   if (fLastTokenStart==0) {//skip any starting delimiters
     delpos=(char*)chars();   
     while (*delpos!='\0' && strchr(fTokenDelimiter, *delpos)!=NULL) 
       delpos++;
     if (*delpos!='\0')
         fLastTokenStart = delpos-chars();
       else { //only delimiters here,no tokens
         GFREE(fTokenDelimiter);
         fLastTokenStart=0;
         return false;
         }
     }
   //now fLastTokenStart is on a non-delimiter char
   //GMessage("String at fLastTokenStart=%d is %s\n", fLastTokenStart, delpos);
   char* token_end=NULL;
   delpos=(char*)strpbrk(chars()+fLastTokenStart,fTokenDelimiter);
   if (delpos==NULL) delpos=(char*)(chars()+length());
   token_end=delpos-1;
   while (*delpos!='\0' && strchr(fTokenDelimiter, *delpos)!=NULL) 
      delpos++; //skip any other delimiters in the set!
   //now we know that delpos is on the beginning of next token
   tlen=(token_end-chars())-fLastTokenStart+1;
   if (tlen==0) {
       GFREE(fTokenDelimiter);
       fLastTokenStart=0;
       return false;
       }
   token.replace_data(tlen);
   ::memcpy(token.chrs(), &chars()[fLastTokenStart], tlen);
   fLastTokenStart=delpos-chars();
   return true;
   }
 //return true;
}

size_t GString::read(FILE* stream, char* delimiter, size_t bufsize) {
//read up to (and including) the given delimiter string
 if (readbuf==NULL) {
    GMALLOC(readbuf, bufsize);
    readbufsize=bufsize;
    }
  else if (bufsize!=readbufsize) {
            GFREE(readbuf);
            if (bufsize>0) {
              GMALLOC(readbuf, bufsize);
              }
            readbufsize=bufsize;
            }
 if (bufsize==0) {
    replace_data(0);
    return 0; //clear the string and free the buffer
    }
 size_t numread;
 size_t acc_len=0; //accumulated length
 int seplen=strlen(delimiter);
 void* p=NULL;
 Data *data = new_data(0);
 do {
   numread=fread(readbuf, 1, bufsize, stream);
   if (numread) {
     p=Gmemscan(readbuf, bufsize, delimiter, seplen);
     if (p!=NULL) {//found the delimiter
           //position the stream after it
           int l = (char*)p-(char*)readbuf;
           fseek(stream, l+seplen-numread, SEEK_CUR);
           numread=l+seplen;
           }
        else {//not found, go back if not eof
           if (numread==bufsize) {
               fseek(stream, -seplen, SEEK_CUR); //check if this works!
               numread-=seplen;
               }
           }
      if (data==&null_data) {
        data=new_data(numread);
        ::memcpy(data->chars, readbuf, numread);
        acc_len+=numread;
        }
       else {
        GREALLOC(data, sizeof(Data)+acc_len+numread);
        memcpy(&data->chars[acc_len], readbuf, numread);
        acc_len+=numread;
        data->length=acc_len;
        data->chars[acc_len]='\0';
        }
      } //if something read
   } while (p==NULL && numread!=0);
  replace_data(data); 
  return acc_len;
}


int GString::asInt(int base /*=10 */) {
 return strtol(text(), NULL, base);
}

double GString::asReal() {
 return strtod(text(), NULL);
}



int GString::peelInt() const {
 if (is_empty()) return 0;
 char buf[24];
 bool started=false;
 int j=0;
 int i;
 for (i=0;i<length();i++) {
  if (started) {
    if (isdigit(my_data->chars[i])) j++; //set coord
                               else break; //finished
    }
   else
    if (isdigit(my_data->chars[i])) {
        j++; started=true;
        }
  }
 if (j>0) {
   strncpy(buf, &my_data->chars[i-j], j);
   buf[j]='\0';
   return strtol(buf, NULL, 10);
   }
 return 0;
}

int GString::peelIntR() const {
 if (is_empty()) return 0;
 char buf[24];
 bool started=false;
 int j=0;
 int i;
 for (i=length()-1;i>=0;i--) {
  if (started) {
    if (isdigit(my_data->chars[i])) j++; //set length
                               else break; //finished
    }
   else
    if (isdigit(my_data->chars[i])) {
      j++; started=true;
      }
  }
 if (j>0) {
   strncpy(buf, &my_data->chars[i+1], j);
   buf[j]='\0';
   return strtol(buf, NULL, 10);
   }
 return 0;
}

GString GString::to(char c) { //return the first part up to first occurence of c
 int i=index(c);
 if (i>=0) return substr(0,i);
      else return (*this);
}
                           //or whole string if c not found
GString GString::from(char c) { //same as to, but starting from the right side
 int i=rindex(c);
 if (i>=0) return substr(i+1);
      else return (*this);
}

int GString::count(char c){
 //return the number of occurences of char c within the string
 int result=0;
 for (int i=0;i<length();i++)
    if (my_data->chars[i]==c) result++;
 return result;
 }

//=========================================

void GString::invalid_args_error(const char *fname) {
    GError("GString:: %s  - invalid arguments\n", fname);
}

//****************************************************************************

void GString::invalid_index_error(const char *fname) {
    GError("GString:: %s  - invalid index\n", fname);
}
//****************************************************************************

