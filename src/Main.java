
import java.math.*;//TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
import java.text.DecimalFormat;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
public class Main {

    public static  Vector solveR(Matrix R, Vector b, int flag)
    {
        double sum = 0;
        Vector x = new Vector(b.getSize());
    if(flag == 0)
        for(int i = 0;i<b.getSize();++i)
        {
            for(int j = 0;j<i;++j)
            {
                sum+= R.getElemnt(i,j)*x.getElement(j);
            }
            sum = b.getElement(i) - sum;
            sum = sum / R.getElemnt(i,i);
            x.setElement(i,sum);
            sum = 0;
        }
    if(flag == 1)
        for(int i = b.getSize()-1;i>=0;--i)
        {
            for(int j = b.getSize()-1;j>i;--j)
            {
                sum+= R.getElemnt(i,j)*x.getElement(j);
            }
            sum = b.getElement(i) - sum;
            sum = sum / R.getElemnt(i,i);
            x.setElement(i,sum);
            sum = 0;
        }
        return x;
    }
    public static Vector yacobi(Matrix A, double error) {
        Matrix D = new Matrix(A.getSize(),1);
        Matrix Q;
        int im,jm;
        Vector eigenvalyes = new Vector(A.getSize());
        double tan, sin, cos, tau,max;
        int k = 0,numberStep=0;
        D.equal(A);
        while (true) {
            ++k;
            im = 0;
            jm = 1;
            for (int i = 0; i < D.getSize(); ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    if(Math.abs(D.getElemnt(i,j))>Math.abs(D.getElemnt(im,jm)))
                    {
                        im = i;
                        jm = j;
                    }

                }
            }
                    if (Math.abs(D.delta()) < error) {System.out.println("Итераций: "+numberStep+"\n");return D.trace();}

                    if (D.getElemnt(im, jm) != 0.0 && im!=jm) {
                        Q = new Matrix(A.getSize(), 1);
                        tau = (D.getElemnt(im, im) - D.getElemnt(jm, jm)) / D.getElemnt(im, jm);
                        tan = 1 / (tau + Math.signum(tau) * Math.sqrt(tau * tau + 1));
                        cos = 1 / (Math.sqrt(tan * tan + 1));
                        sin = tan * cos;
                        Q.setElement(im, im, cos);
                        Q.setElement(jm, jm, cos);
                        Q.setElement(im, jm, sin);
                        Q.setElement(jm, im, -sin);
                        D.equal(Q.multMatrix(D.multMatrix(Q.transposition())));
                        numberStep++;

                    }

                }
            }

    public static Vector QR(Matrix A, double error)
    {
        Matrix B1 = new Matrix(A.getSize(),0);
        Matrix B = new Matrix(A.getSize());
        Vector b = new Vector(A.getSize());
        B1.equal(A);
        Matrix prom = new Matrix(B1.getSize());
        Matrix Hes = new Matrix(A.getSize());
        Vector w = new Vector(B1.getSize()-1);
        Vector a = new Vector(B1.getSize()-1);
        Vector u = new Vector(a.getSize(), 0);

            for(int j = 0;j<a.getSize();++j)
            {
                a.setElement(j,B1.getElemnt(j+1,1));
            }
            System.out.println(a);
            u.setElement(0,a.norm());
        System.out.println(u);
        Hes.equal(B1.sumMatrix(new Matrix(a.getSize(),1), a.uut(a,a).multOnScal(2/(a.utu(a,a)))));
        for(int i = 0;i<prom.getSize();++i)
        {
            for (int j = 0;j<prom.getSize();++j)
            {
                if(i == 0 || j ==0)
                {
                    prom.setElement(i,j,0);
                }
                else
                {
                    prom.setElement(i,j,Hes.getElemnt(i-1,j-1));
                }
                Hes.setElement(0,0,B1.getElemnt(0,0));
            }
        }
        int numberStep = 0;
        while (true)
        {
            B.equal(rollingMethod11(B1,b));
            if(B1.sumdown()<error)
            {
                System.out.println("Итераций: "+numberStep+"\n");
                return B.trace();
            }
            B1.equal(B);
            ++numberStep;
        }
    }
    public static Vector conjugateGradient(Matrix A,Vector b,double epsilon)
    {
        int k = 0;
        double lamda = Math.pow(QR(A, epsilon).chebNorm1(), -1), step, v;
        Vector x = new Vector(A.getSize(), 0);
        Vector discrepancy = b.sumVector(A.multOnVector(x).multOnScal(-1.0));
        Vector discrepancy1 = new Vector(b.getSize());
        discrepancy1.equal(discrepancy);
        Vector s = new Vector(b.getSize());
        Vector g = new Vector(b.getSize());
        s.equal(discrepancy);
        try {
            File fac = new File("C:\\work\\Git_Projects\\VCH4\\graphic.txt");
            FileWriter wr = new FileWriter(fac);
            while (lamda * (A.multOnVector(x).sumVector(b.multOnScal(-1))).norm() > epsilon) {
                g.equal(A.multOnVector(s));
                step = discrepancy.scalMult(discrepancy) / s.scalMult(g);
                x.equal(x.sumVector(s.multOnScal(step)));
                discrepancy.equal(discrepancy.sumVector(g.multOnScal(step * (-1))));
                v = discrepancy.scalMult(discrepancy) / discrepancy1.scalMult(discrepancy1);
                s.equal(discrepancy.sumVector(s.multOnScal(v)));
                discrepancy1.equal(discrepancy);
                k++;
                wr.write(Integer.toString(k)+","+Double.toString(Math.log10(A.multOnVector(x).sumVector(b.multOnScal(-1.0)).norm()))+
                        System.getProperty( "line.separator" ));
            }
            wr.close();
            System.out.println("Итераций: "+k+"\n");
            return x;
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }
        return x;
    }
    public static Vector fallDawn(Matrix A,Vector b,double epsilon)
    {
        Matrix D = new Matrix(A.getSize());
        D.equal(A);
        double step,lamda = Math.pow(QR(A,epsilon).chebNorm1(),-1);
        Vector x = new Vector(b.getSize(),1);
        Vector discrepancy = new Vector(b.getSize(),0);
        int numberStep =0;

        //(0.5*A.multOnVector(x).scalMult(x) - b.scalMult(x))<=epsilon
        try {
            File fac = new File("C:\\work\\Git_Projects\\VCH4\\graphic.txt");
            FileWriter wr = new FileWriter(fac);

        while (lamda*(A.multOnVector(x).sumVector(b.multOnScal(-1))).norm()>epsilon && numberStep<100000)
        {
            discrepancy.equal(A.multOnVector(x).sumVector(b.multOnScal(-1.0)));
            step = ((Math.pow(discrepancy.norm(),2))/(A.multOnVector(discrepancy).scalMult(discrepancy)));
            x.equal(x.sumVector(discrepancy.multOnScal(-step)));
            numberStep++;
            wr.write(Integer.toString(numberStep)+","+Double.toString(Math.log10(A.multOnVector(x).sumVector(b.multOnScal(-1.0)).norm()))+
                    System.getProperty( "line.separator" ));
        }
            wr.close();
            System.out.println("Итераций: "+numberStep+"\n");
            return x;
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }
        return x;
    }
    public static Matrix rollingMethod11(Matrix A,Vector b)
    {
        Matrix Q = new Matrix(A.getSize(),1);
        Matrix R = new Matrix(A.getSize(),1);
        Matrix roll;
        R.equal(A);
        double cosa,sina;
        for(int i = 0;i<A.getSize();++i)
        {
            for (int j = i+1;j<A.getSize();++j)
            {
                roll = new Matrix(A.getSize(), 1);
                cosa = R.getElemnt(i,i)/(Math.sqrt(R.getElemnt(i,i)*R.getElemnt(i,i) +
                        R.getElemnt(j,i)*R.getElemnt(j,i)));
                sina = (-R.getElemnt(j,i))/(Math.sqrt((R.getElemnt(i,i)*R.getElemnt(i,i) +
                        R.getElemnt(j,i)*R.getElemnt(j,i))));
                roll.setElement(i,i,cosa);
                roll.setElement(j,i,sina);
                roll.setElement(i,j,-sina);
                roll.setElement(j,j,cosa);
                R.equal(roll.multMatrix(R));
                Q.equal(roll.multMatrix(Q));

            }
        }
        return R.multMatrix(Q.transposition());
    }
    public static Vector rollingMethod(Matrix A,Vector b)
    {
        Matrix Q = new Matrix(A.getSize(),1);
        Matrix R = new Matrix(A.getSize(),1);
        Matrix roll;
        R.equal(A);
        double cosa,sina;
        for(int i = 0;i<A.getSize();++i)
        {
         for (int j = i+1;j<A.getSize();++j)
         {
             roll = new Matrix(A.getSize(), 1);
             cosa = R.getElemnt(i,i)/(Math.sqrt(R.getElemnt(i,i)*R.getElemnt(i,i) +
                                            R.getElemnt(j,i)*R.getElemnt(j,i)));
             sina = (-R.getElemnt(j,i))/(Math.sqrt((R.getElemnt(i,i)*R.getElemnt(i,i) +
                                            R.getElemnt(j,i)*R.getElemnt(j,i))));
             roll.setElement(i,i,cosa);
             roll.setElement(j,i,sina);
             roll.setElement(i,j,-sina);
             roll.setElement(j,j,cosa);
             R.equal(roll.multMatrix(R));
             Q.equal(roll.multMatrix(Q));

         }
        }
        System.out.println(R);
        System.out.println(Q);
        System.out.println("//////////////////////");
        return solveR(R,Q.multOnVector(b),1);
    }
    public static double rollingMethodmin(Matrix A)
    {
        Matrix R = new Matrix(A.getSize(),1);
        Matrix roll;
        R.equal(A);
        double cosa,sina,min = 0;
        for(int i = 0;i<A.getSize();++i)
        {
            for (int j = i+1;j<A.getSize();++j)
            {
                roll = new Matrix(A.getSize(), 1);
                cosa = R.getElemnt(i,i)/(Math.sqrt(R.getElemnt(i,i)*R.getElemnt(i,i) +
                        R.getElemnt(j,i)*R.getElemnt(j,i)));
                sina = (-R.getElemnt(j,i))/(Math.sqrt((R.getElemnt(i,i)*R.getElemnt(i,i) +
                        R.getElemnt(j,i)*R.getElemnt(j,i))));
                roll.setElement(i,i,cosa);
                roll.setElement(j,i,sina);
                roll.setElement(i,j,-sina);
                roll.setElement(j,j,cosa);
                R.equal(roll.multMatrix(R));

            }
        }
        min = R.getElemnt(0,0);
        for(int i = 0;i< R.getSize();++i)
        {
            if(min< R.getElemnt(i,i))
            min = R.getElemnt(i,i);
        }
        return Math.abs(min);
    }
    public static Vector holess(Matrix A, Vector b)
    {
        double sum = 0;
        Vector y;
        Vector x;
        Matrix hol = new Matrix(A.getSize(),0);

        for(int j = 0;j<hol.getSize();j++)
        {
            sum = 0;
            for(int k = 0;k<j;++k)
            {
                sum += hol.getElemnt(j,k) * hol.getElemnt(j,k);
            }
            sum = A.getElemnt(j,j) - sum;
            if(sum<0)
            {
                System.out.printf("Попытка извлечь корень из отрицательного %.16e\n", sum);
                break;
            }
            sum = Math.sqrt((sum));
            hol.setElement(j,j,sum);
            sum = 0;
            for(int i = j+1;i< A.getSize();i++)
            {
                for(int k = 0;k<j;++k)
                {
                    sum += hol.getElemnt(i,k) * hol.getElemnt(j,k);
                }
                sum = A.getElemnt(i,j) - sum;
                sum = sum / hol.getElemnt(j,j);
                hol.setElement(i,j,sum);
                sum = 0;

            }
        }
        y = solveR(hol,b,0);
        x = solveR(hol.transposition(),y,1);
        return x;
    }
    public static void main(String[] args) {
        int n =15, k = 13;
        DecimalFormat df = new DecimalFormat("0.000000000000000E0");
        double iterations = 0.0000001;
        double time,timeall=0;
        double[] vector_massive = new double[n];
        for(int i = 0;i<n;++i)
        {
            vector_massive[i] = Math.pow(-1.0,i);
        }
        Matrix Gilbert = new Matrix(n,2);
        Matrix Gilbert1 = new Matrix(n,3);
        Matrix A = new Matrix(n,3);
        Vector answers = new Vector(vector_massive,n);
        System.out.println(A);
//        for(int i = 0;i<100;++i)
//        {
//            A.equal(rollingMethod11(A,answers));
//        }
//
//        System.out.println(Gilbert);
//        System.out.println(answers);
//        //System.out.println(Gilbert.determinant(Gilbert1,n));
//        time = System.currentTimeMillis();
//        System.out.println("Holess\n");
//        System.out.println("Answer = " + holess(Gilbert1,answers));
//        System.out.println("Left part = " + Gilbert.multOnVector(holess(Gilbert1,answers)));
//        System.out.println("Run time = " + (System.currentTimeMillis() - time)/2 + "\n");
//        System.out.println("Nevyazka = " + Gilbert.multOnVector(holess(Gilbert1,answers)).sumVector(answers.multOnScal(-1.0)));
//        System.out.println("Norm of nevyazka = " + (Gilbert.multOnVector(holess(Gilbert,answers)).sumVector(answers.multOnScal(-1.0))).norm());
//        System.out.println("\n\n _____________________________ \n\n");
//        System.out.println("Rolling method\n");
//        time = System.currentTimeMillis();
//        System.out.println("Answer = " + rollingMethod(Gilbert1,answers));
//        System.out.println("Left part = " + Gilbert.multOnVector(rollingMethod(Gilbert1,answers)));
//        System.out.println("Run time = " + (System.currentTimeMillis() - time)/2);
//        System.out.println("Nevyazka = " + Gilbert.multOnVector(rollingMethod(Gilbert1,answers)).sumVector(answers.multOnScal(-1.0)));
//        System.out.printf("Norm of nevyazka = %.16e\n", (Gilbert.multOnVector(rollingMethod(Gilbert1,answers)).sumVector(answers.multOnScal(-1.0))).norm());
//        System.out.printf("Разница = %.16e\n", (Gilbert.multOnVector(rollingMethod(Gilbert1,answers)).sumVector(answers.multOnScal(-1.0)).norm()
//                - Gilbert.multOnVector(holess(Gilbert,answers)).sumVector(answers.multOnScal(-1.0)).norm()));
////        System.out.println("Falling Dawn method\n");
//        System.out.println("Answer = " + fallDawn(Gilbert1,answers,iterations));
//        System.out.println("Left part = " + Gilbert.multOnVector(fallDawn(Gilbert1,answers,iterations)));
//        System.out.println("Run time = " + (System.currentTimeMillis() - time));
//        System.out.println("Nevyazka = " + Gilbert.multOnVector(fallDawn(Gilbert1,answers,iterations)).sumVector(answers.multOnScal(-1.0)));
//        System.out.println("Norm of nevyazka = " + (Gilbert.multOnVector(fallDawn(Gilbert1,answers,iterations)).sumVector(answers.multOnScal(-1.0))).norm());
//        System.out.println("\n Разница = " + (Gilbert.multOnVector(fallDawn(Gilbert1,answers,iterations)).sumVector(answers.multOnScal(-1.0)).norm()
//                - Gilbert.multOnVector(holess(Gilbert,answers)).sumVector(answers.multOnScal(-1.0)).norm()));
      //  time = 0;
//        System.currentTimeMillis();
//        for(int i = 0;i<k;++i)
//        {
//            time = System.currentTimeMillis();
//             Vector c = holess(Gilbert, answers);
//             time = System.currentTimeMillis() - time;
//             timeall += time;
//        }
//        System.out.printf("Run time Holess = %.16e\n", timeall/k);
//        timeall =0;
//        for(int i = 0;i<k;++i)
//        {
//            time = System.currentTimeMillis();
//            Vector c = rollingMethod(Gilbert, answers);
//            time = System.currentTimeMillis() - time;
//            timeall += time;
//        }
//        System.out.printf("Run time Roll = %.16e\n", timeall/k);


//        System.out.println(A);
//        System.out.println(yacobi(A,iterations));
//        System.out.println(QR(A,iterations));
//        System.out.println(yacobi(A,iterations).sumVector(QR(A,iterations).multOnScal(-1)).norm());
//        System.out.println("***************************************************************");
//        System.out.println(holess(A,answers));
        System.out.println(A);
        System.out.println("Собсьвенные значения вычесленные методом Якоби :\n");
        System.out.println(yacobi(A,iterations));
        System.out.println("Собсьвенные значения вычесленные методом QR :\n");
        System.out.println(QR(A,iterations));
        Vector eiginevalues = QR(A,iterations);

//        System.out.println("метод наскорейшего градиентного спуска: \n");
//        System.out.println(fallDawn(A,answers,iterations));
//        System.out.println("Невязка: \n");
//        System.out.println("'''"+A.multOnVector(fallDawn(A,answers,iterations)).sumVector(answers.multOnScal(-1.0)
//        ).norm());
//        System.out.println("Метод сопряженных градиентов: \n");
//        System.out.println(conjugateGradient(A,answers,iterations));
//        System.out.println("Невязка: \n");
//        System.out.println("'''"+A.multOnVector(conjugateGradient(A,answers,iterations)).sumVector(answers.multOnScal(-1.0)
//        ).norm());
//        System.out.println("Холесский: \n");
//        System.out.println(holess(A,answers));
                //        System.out.println(A);

//        System.out.println(Gilbert1.sumMatrix(A,
//                new Matrix(eiginevalues.getSize(),1)));
//        System.out.println(new Matrix(eiginevalues.getSize(),1));
//        System.out.println(A);
//        //System.out.println(fallDawn(A,answers,iterations));
//        System.out.println(conjugateGradient(A,answers,iterations));
//        //System.out.println(fallDawn(A,answers,iterations));
//        //System.out.println(rollingMethod(A.sumMatrix(new Matrix(n,3),new Matrix(n,1).multOnScal(-1*eiginevalues.getElement(1))),new Vector(n,0)));
//        System.out.println((A.sumMatrix(new Matrix(n,3),new Matrix(n,1).multOnScal(-1*eiginevalues.getElement(1)))));



    }
}