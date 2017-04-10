

#include "segment.h"
#include "reflective.h"

reflective::reflective(double x1, double y1, double x2, double y2, double sp){
    point1.x = x1; point1.y = y1; point2.x = x2; point2.y = y2; spec = sp;
    length = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    if (length!=0)
    {
        normal.x = -(y2-y1)/length;
        normal.y = (x2-x1)/length;
    }
    else{
        normal.x = 0;
        normal.y = 0;
    }
}
