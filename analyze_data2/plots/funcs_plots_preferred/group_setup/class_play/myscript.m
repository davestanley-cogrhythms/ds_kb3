

clear 
a(1:10) = MySuperClass
b(1:10) = MySubClass
%b(1) = a(1)

b(1).inheritObj(a(1));

a(2).inheritObj(b(2));
