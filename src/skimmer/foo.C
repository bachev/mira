#include <iostream>

class FOO
{
private:
  void foofunc(unsigned short * bar);
public:
  void dofoo();
};

void FOO::foofunc(unsigned short * bar)
{
  cout << "1";
}

void FOO::dofoo()
{
  unsigned short * arr = new unsigned short[10000];
  foofunc(arr);
}

int main()
{
  FOO test;
  test.dofoo();
}
