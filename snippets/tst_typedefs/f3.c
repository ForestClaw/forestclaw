struct s;

typedef void (*f_t)(struct s);

typedef struct s
{
    f_t g;
} s_t;


int main()
{
    return 0;
}
