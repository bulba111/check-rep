#include <iostream>
using namespace std;

class myclass
{
public:
    myclass(int size)
    {
        this-> Size = size;
        data = new int(size);
        for (int i = 0; i<size; i++) {
            data[i] = i;
        }
        
        cout << "constructor " << this << endl;
    }
    
    myclass(const myclass& other)
    {
        this->Size = other.Size;
        
        this->data = new int[other.Size];
        for (int i = 0; i<other.Size; i++)
        {
            this->data[i] = other.data[i];
        }
        
        
        cout << "copy constructor " << this << endl;
    }
    
    ~myclass()
    {
        delete[] data;
        cout << "destructor " << this << endl;
    }

private:
    int* data;
    int Size;
};



int main()
{
    myclass a(5);
    myclass b(a);
    return 0;
}
