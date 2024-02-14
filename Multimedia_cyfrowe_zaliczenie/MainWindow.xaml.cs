using OxyPlot.Series;
//using OxyPlot.Wpf;
using OxyPlot;
using System;
using System.Linq;
using System.Windows;
using System.Windows.Input;
using System.Data;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;
using OxyPlot.Axes;

namespace Multimedia_cyfrowe_zaliczenie
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    /// Application made by Mikolaj Luczak
    public partial class MainWindow : System.Windows.Window
    {
        private PlotModel plotModel;

        //range
        private bool isMultipliedByPi = true;
        private double minValue = 0;
        private double maxValue = 2 * Math.PI;
        
        //number of points
        private int numberOfPoints = 128;

        //first signal parameters
        private double firstSignalFrequency = 1;
        private double firstSignalAmplitude = 1;
        private double firstSignalPhase = 0;
        private string firstSignalType = "SINUS";

        //second signal parameters
        private double secondSignalFrequency = 1;
        private double secondSignalAmplitude = 1;
        private double secondSignalPhase = 0;
        private string secondSignalType = "SINUS";

        //second signal parameters
        private double thirdSignalFrequency = 1;
        private double thirdSignalAmplitude = 1;
        private double thirdSignalPhase = 0;
        private string thirdSignalType = "SINUS";

        //additional parameters
        private bool isNoiseAdded = false;
        private bool isFFTCalculated = false;
        private bool isWindowSelected = false;
        private string windowType = "BLACKMAN";

        private int numberOfPointsFFT = 16;

        //signals tab for fourier transform
        private double[] resultSignal;
        private double[] fftResultSignal;

        public MainWindow()
        {
            InitializeComponent();
        }

        //function calculating window
        private double AddWindow(double n)
        {
            double N = numberOfPoints;
            double result = 0;
            switch (windowType)
            {
                case "BLACKMAN":
                    result = 0.42 * 0.5 * Math.Cos((2 * Math.PI * n) / (N - 1)) - 0.08 * Math.Cos(4 * Math.PI * n / (N - 1));
                    break;
                case "HANNIMG":
                    result = 0.5 * (1 - Math.Cos((2 * Math.PI * n) / N));
                    break;
                case "HAMMING":
                    result = 0.54 - 0.46 * Math.Cos((2 * Math.PI * n) / N);
                    break;
                case "BARLET":
                    result = 1 - Math.Abs(n - (N - 1) / 2) / ((N - 1) / 2);
                    break;
                default:
                    result = n;
                    break;
            }

            return result;
        }
        
        //function calculating value of signal in point
        private double calculateY (double x, double frequency, double amplitude, double phase,string signalType)
        {
            double y = 0;
            switch(signalType)
            {
                case "SINUS":
                    y = amplitude * Math.Sin(frequency * (x + phase)); 
                    break;
                case "SQUARE":
                    double sinx = amplitude * Math.Sin(frequency * (x + phase));
                    y = (sinx >= 0) ? amplitude : -1 * amplitude;
                    break;
                case "SAWTOOTH":
                    y = amplitude * ((frequency * (x + phase)) % 1); 
                    break;
                default:
                    y = 0; 
                    break;
            }

            if (isNoiseAdded)
            {
                Random random = new Random();
                y += ((random.NextDouble() * 0.06) - 0.03);//* amplitude; //??
            }

            return y;
        }

        //function that calculate and draw signal
        private void ShowPlot()
        {
            // Initialization of plot model 
            plotModel = new PlotModel();
            plotModel.Title = "Wykres";

            // Adding line series to the model
            LineSeries signalSeries = new LineSeries();
            //signalSeries.Title = "Result";

            //calculation every point of signal based on chose parameters
            for (int i = 0; i < numberOfPoints; i++)
            {
                double x = minValue + (i * (maxValue - minValue)) / (numberOfPoints - 1); //zakres zmienny
                //double x = i * (2 * Math.PI) / (numberOfPoints - 1); //zakres staly od 0 do 2pi
                double y1 = calculateY(x, firstSignalFrequency, firstSignalAmplitude, firstSignalPhase,firstSignalType);
                double y2 = calculateY(x, secondSignalFrequency, secondSignalAmplitude, secondSignalPhase,secondSignalType);
                double y3 = calculateY(x, thirdSignalFrequency, thirdSignalAmplitude, thirdSignalPhase, thirdSignalType);
                double y = y1 + y2 + y3;
                if (isWindowSelected && !isFFTCalculated)
                    y = y * AddWindow(x);
                signalSeries.Points.Add(new DataPoint(x, y));
            }

            resultSignal = signalSeries.Points.Select(p => p.Y).ToArray();
            fftResultSignal = calculateFourierTransform(resultSignal);
            
            //calculating fourier transform - if we selected this option
            //LineSeries fftSeries = new LineSeries();
            ColumnSeries fftSeries = new ColumnSeries();
            for (int i = 0; i < numberOfPointsFFT; i++)
            {
                int x = i;
                double y = fftResultSignal[i];
                if (isWindowSelected)
                    y = y * AddWindow(x);
                y = Math.Abs(fftResultSignal[i]); // / numberOfPointsFFT;
                //fftSeries.Points.Add(new DataPoint(x, y));
                fftSeries.Items.Add(new ColumnItem(y, x - 1)); //new ColumnItem(y);
            }

            // Clearing old plot and adding new one
            plotModel.Series.Clear();
            if (isFFTCalculated)
            {
                plotModel.Series.Add(fftSeries);

                var categoryAxis = new CategoryAxis { Position = AxisPosition.Bottom };
                categoryAxis.IntervalLength = 1;
                categoryAxis.GapWidth = 0.1;
                int step = (numberOfPointsFFT / 100).Round(0) * 10; // step is 1 for 10 point 10 for 100 points etc.
                if (step < 1)
                    step = 1;
                categoryAxis.MajorStep = step;
                var linearAxis = new LinearAxis { Position = AxisPosition.Right };

                plotModel.Axes.Add(categoryAxis);
                plotModel.Axes.Add(linearAxis);
            }
            else
            {
                plotModel.Series.Add(signalSeries);
            }

            // Assinging model to PlotView from GUI
            plotView.Model = plotModel;
        }
        
        // function that calculate fourier transform
        // it takes and return complex numbers so we have to take only real part of the number
        private double[] calculateFourierTransform(double[] signal)
        {
            Complex[] complexSignal = Array.ConvertAll(signal, x => new Complex((float)x, 0.0f));

            Fourier.Forward(complexSignal, FourierOptions.NoScaling);

            double[] result = new double[complexSignal.Length];

            for (int i = 0; i < complexSignal.Length; i++)
            {
                result[i] = complexSignal[i].Real;
            }

            return result;
        }
        
        // we need to get values from UI and assign it to parameters
        private void drawButton_Click(object sender, RoutedEventArgs e)
        {
            isMultipliedByPi = multiplyByPiCheckBox.IsChecked.Value;

            minValue = isMultipliedByPi ? double.Parse(minimalValueTextInput.Text) * Math.PI   : double.Parse(minimalValueTextInput.Text);
            maxValue = isMultipliedByPi ? double.Parse(maximalValueTextInput.Text) * Math.PI : double.Parse(maximalValueTextInput.Text);
            numberOfPoints = int.Parse(numberOfPointsTextInput.Text);

            firstSignalFrequency = double.Parse(firstSignalFrequencyTextInput.Text);
            firstSignalAmplitude = double.Parse(firstSignalAmplitudeTextInput.Text);
            firstSignalPhase = double.Parse(firstSignalPhaseTextInput.Text);
            if (firstSignalTypeDropDown.SelectedItem != null)
                firstSignalType = firstSignalTypeDropDown.SelectedItem.ToString().Split(' ')[1];

            secondSignalFrequency = double.Parse(secondSignalFrequencyTextInput.Text);
            secondSignalAmplitude = double.Parse(secondSignalAmplitudeTextInput.Text);
            secondSignalPhase = double.Parse(secondSignalPhaseTextInput.Text);
            if (secondSignalTypeDropDown.SelectedItem != null)
                secondSignalType = secondSignalTypeDropDown.SelectedItem.ToString().Split(' ')[1];

            thirdSignalFrequency = double.Parse(thirdSignalFrequencyTextInput.Text);
            thirdSignalAmplitude = double.Parse(thirdSignalAmplitudeTextInput.Text);
            thirdSignalPhase = double.Parse(thirdSignalPhaseTextInput.Text);
            if (thirdSignalTypeDropDown.SelectedItem != null)
                thirdSignalType = thirdSignalTypeDropDown.SelectedItem.ToString().Split(' ')[1];

            numberOfPointsFFT = int.Parse(numberOfPointsFFTTextInput.Text);
            if (numberOfPointsFFT > numberOfPoints)
                numberOfPointsFFT = numberOfPoints;

            isNoiseAdded = noiseCheckBox.IsChecked.Value;
            isFFTCalculated = fourierCheckBox.IsChecked.Value;
            isWindowSelected = windowCheckBox.IsChecked.Value;
            if (firstSignalTypeDropDown.SelectedItem != null)
                windowType = windowTypeDropDown.SelectedItem.ToString().Split(' ')[1];
            
            ShowPlot();

            //MessageBox.Show($"ilosc punktow {numberOfPoints} \n czest1 {firstSignalFrequency} \n amp1 {firstSignalAmplitude} \n czest2 {secondSignalFrequency} \n amp2 {secondSignalAmplitude} \n punkty fft {numberOfPointsFFT}");
            //MessageBox.Show($"1sig {firstSignalType}, 2sig {secondSignalType}");
            //MessageBox.Show($"{isNoiseAdded}");
            //MessageBox.Show($" {windowType}");
            //MessageBox.Show($" {isMultipliedByPi} {minValue} {maxValue}");
        }

        private void TextBox_IntegerOnly(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !IsNumeric(e.Text);
        }

        private bool IsNumeric(string text)
        {
            return int.TryParse(text, out _);
        }

        private void TextBox_DoubleOnly(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !IsDouble(e.Text);
        }

        private bool IsDouble(string text)
        {
            return double.TryParse(text, out _) || text == ",";
        }
    }
}
