// functions for hashing
#include <tuple>

namespace std
{

template<typename T>
struct hash<std::tuple<T, T>>
{
    size_t operator()(const std::tuple<T, T>& arg) const
    {
        size_t h1 = std::hash<T>{}(std::get<0>(arg));
        size_t h2 = std::hash<T>{}(std::get<1>(arg));
        return h1 ^ (h2 << 1);
    }
};

template<typename T>
struct hash<std::tuple<T, T, T>>
{
    size_t operator()(const std::tuple<T, T, T>& arg) const
    {
        size_t h1 = std::hash<T>{}(std::get<0>(arg));
        size_t h2 = std::hash<T>{}(std::get<1>(arg));
        size_t h3 = std::hash<T>{}(std::get<2>(arg));
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};
}
