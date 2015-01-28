clear;
x = 3; 
y = 4;
fun = @(x, y, z) x * 100 + 10 * y + z;
for z = 1:4
    fun(1,2,z)
end