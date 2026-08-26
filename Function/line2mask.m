function BW = line2mask(x, y, M, N)

BW = zeros(M, N);

p1Outside = x(1) < 1 || x(1) > N || y(1) < 1 || y(1) > M;
p2Outside = x(2) < 1 || x(2) > N || y(2) < 1 || y(2) > M;
if p1Outside && p2Outside, return; end

[r0, c0, r1, c1] = clipPoints(x(:)', y(:)', M, N);

BW(r0, c0) = 1;
BW(r1, c1) = 1;
if abs(r1 - r0) <= abs(c1 - c0)
   if c1 < c0
      k = r1; r1 = r0; r0 = k;
      k = c1; c1 = c0; c0 = k;
   end
   if (r1 >= r0) & (c1 >= c0)
      dy = c1-c0; dx = r1-r0;
      p = 2*dx; n = 2*dy - 2*dx; tn = dy;
      while (c0 < c1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; r0 = r0 + 1;
         end
         c0 = c0 + 1; BW(r0, c0) = 1;
      end
   else
      dy = c1 - c0; dx = r1 - r0;
      p = -2*dx; n = 2*dy + 2*dx; tn = dy;
      while (c0 <= c1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; r0 = r0 - 1;
         end
         c0 = c0 + 1; BW(r0, c0) = 1;
      end
   end
else if r1 < r0
      k = r1; r1 = r0; r0 = k;
      k = c1; c1 = c0; c0 = k;
   end
   if (r1 >= r0) & (c1 >= c0)
      dy = c1 - c0; dx = r1 - r0;
      p = 2*dy; n = 2*dx-2*dy; tn = dx;
      while (r0 < r1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; c0 = c0 + 1;
         end
         r0 = r0 + 1; BW(r0, c0) = 1;
      end
   else
      dy = c1 - c0; dx = r1 - r0;
      p = -2*dy; n = 2*dy + 2*dx; tn = dx;
      while (r0 < r1)
         if tn >= 0
            tn = tn - p;
         else
            tn = tn + n; c0 = c0 - 1;
         end
         r0 = r0 + 1; BW(r0, c0) = 1;
      end
   end
end
BW = BW(1:M, 1:N);
end
