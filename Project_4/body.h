#ifndef BODY_H
#define BODY_H


class body
{
public:

    double mass;
            double position[4];
            double velocity[4];
            body();
            body(double, double, double, double, double, double, double);
            /*body(body, double, double, double, double, double, double, double);*/
};

#endif // BODY_H
