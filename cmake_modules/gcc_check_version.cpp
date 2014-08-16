#if __GNUC__  < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 6)
        #error Minimum GCC version supported is 4.6.0.
#endif

int main()
{
    return 0;
}
