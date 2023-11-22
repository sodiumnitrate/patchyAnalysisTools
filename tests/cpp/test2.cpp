#include "vec3.hpp"

int main()
{
    Vec3 vec1(0,0,0);
    Vec3 vec2(2,2,2);
    Vec3 vec3(3,3,5);

    Vec3 disp;

    disp = vec3.diff(vec2);

    //return (disp.get_x() == 1 && disp.get_y() == 1 && disp.get_z() == 3);
    return !vec2.is_normal();
}