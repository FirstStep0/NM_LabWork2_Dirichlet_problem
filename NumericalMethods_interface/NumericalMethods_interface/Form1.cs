using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using Newtonsoft.Json;

namespace NumericalMethods_interface
{
    public partial class Form1 : Form
    {
        Random rand = new Random();
        DataTable values;

        global_data data = new global_data();

        string path_prog = "NumericalMethods_Lab2.exe";
        string path_input = "input.txt";
        string path_output = "result.txt";
        string path_graph = "3d_graph.py";
        public Form1()
        {
            InitializeComponent();
        }
        private bool parseParams(global_data data) {
            data.numberVariant = 1;
           // data.numberMethod = 1;
            data.a = 1;
            data.b = 2;
            data.c = 2;
            data.d = 3;
            data.area = "1";
            data.area_x = data.area_y = 1;
            int nmax1, nmax2;
            double eps1, eps2;
            double param;
            if (!int.TryParse(textBox_Nmax1.Text, out nmax1)) MessageBox.Show("Incorrect N_max1");
            else if (!int.TryParse(textBox_Nmax2.Text, out nmax2)) MessageBox.Show("Incorrect N_max2");
            else if (!double.TryParse(textBox_eps1.Text, out eps1)) MessageBox.Show("Incorrect eps1");
            else if (!double.TryParse(textBox_eps2.Text, out eps2)) MessageBox.Show("Incorrect eps2");
            else if (!int.TryParse(textBox_n.Text, out data.n)) MessageBox.Show("Incorrect n");
            else if (!int.TryParse(textBox_m.Text, out data.m)) MessageBox.Show("Incorrect m");
            else if (comboBox2.SelectedIndex < 0) MessageBox.Show("Incorrect task");
            else if (!double.TryParse(textBox_param.Text, out param)) MessageBox.Show("Incorrect param");
            else
            {
                data.numberMethod = comboBox1.SelectedIndex + 1;
                data.nmax.Add(nmax1);
                data.nmax.Add(nmax2);
                data.eps.Add(eps1);
                data.eps.Add(eps2);
                data.param.Add(param);
                data.test = (comboBox2.SelectedIndex == 0);
                return true;
            }
            return false;
        }
        private void button1_Click(object sender, EventArgs e)
        {
            try
            {
                data = new global_data();
                if (parseParams(data)) {
                    string json = JsonConvert.SerializeObject(data);
                    FileStream fs = File.Open(path_input, FileMode.Create, FileAccess.Write);
                    fs.SetLength(0);
                    StreamWriter writer = new StreamWriter(fs);
                    writer.Write(json);
                    writer.Close();

                    if (backgroundWorker1.IsBusy) throw new Exception("backgroundWorker is busy");
                    else backgroundWorker1.RunWorkerAsync();
                }
            }
            catch (Exception exc)
            {
                MessageBox.Show(exc.Message);
            }
        }
        private void chart1_MouseClick(object sender, MouseEventArgs e)
        {
            //chart1.Refresh();
        }

        private int Clamp(int x, int min, int max) {
            if (x < min) x = min;
            if (x > max) x = max;
            return x;
        }

        private void backgroundWorker1_DoWork(object sender, DoWorkEventArgs e)
        {
            backgroundWorker1.ReportProgress(0);
            ProcessStartInfo info = new ProcessStartInfo();
            info.FileName = path_prog;
            info.Arguments = path_input + " " + path_output;
            info.RedirectStandardError = true;
            info.RedirectStandardOutput = true;
            info.CreateNoWindow = true;
            info.UseShellExecute = false;
            StreamReader sr;
            Task<string> task = null;


            Process p = Process.Start(info);
            sr = p.StandardOutput;
            while ((!p.HasExited))
            {
                if (backgroundWorker1.CancellationPending) { p.Kill(); e.Cancel = true; return; };
                if (task == null) task = p.StandardOutput.ReadLineAsync();
                if (task != null && task.IsCompleted) {
                    int stage;
                    if (int.TryParse(task.Result, out stage)) {
                        backgroundWorker1.ReportProgress(Clamp((int)(stage), 0, 100));
                    }
                    task = null;
                }
            }
            p.WaitForExit();
            if (p.ExitCode != 0) throw new Exception(p.StandardError.ReadLine());

            backgroundWorker1.ReportProgress(100);
        }

        private void print_ref_for_test_task(answer ans) {
            textBox_ref.Text = "";
            textBox_ref.Text += String.Format("Справка для тестовой задачи\r\n");
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("«Для решения тестовой задачи использованы сетка\r\n");
            textBox_ref.Text += String.Format("с числом разбиений по x n = {0}\tи числом разбиений по y m = {1},\r\n", ans.n, ans.m);
            textBox_ref.Text += String.Format("Использован {0} с параметром ω = {1}\r\n", "метод верхней релаксации", ans.param[0]);
            textBox_ref.Text += String.Format("Применены критерии остановки:\r\n");
            textBox_ref.Text += String.Format(" по точности εмет = {0:e2}\r\n", ans.eps[0]);
            textBox_ref.Text += String.Format(" по числу итераций Nmax = {0}\r\n", ans.nmax[0]);
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("На решение схемы (СЛАУ) затрачено итераций N = {0} \r\n", ans.N[0]);
            textBox_ref.Text += String.Format("и достигнута точность итерационного метода ε(N) = {0:e3}\r\n", ans.acc[0]);
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("Схема (СЛАУ) решена с невязкой ||R(N)|| = {0:e3} по {1}\r\n", ans.R[0], "норме max");
            textBox_ref.Text += String.Format("для невязки СЛАУ использована {0}\r\n", "норма max");
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("Тестовая задача должна быть решена с погрешностью не более ε = 0.5*10^–6;\r\n");
            textBox_ref.Text += String.Format("задача решена с погрешностью ε1 = {0:e3}\r\n", ans.err);
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("Максимальное отклонение точного и численного решений наблюдается в узле x = {0:f3}; y = {1:f3}\r\n", ans.x_err, ans.y_err);
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("В качестве начального приближения использовано {0}\r\n", "интерполяция по x");
            textBox_ref.Text += String.Format("\r\nВремя работы: {0} сек", ans.time);
        }

        private void print_ref_for_basic_task(answer ans) {
            textBox_ref.Text = "";
            textBox_ref.Text += String.Format("Справка для основной задачи\r\n");
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("«Для решения основной задачи использована сетка с числом разбиений по x\r\n");
            textBox_ref.Text += String.Format("n = {0} и числом разбиений по y m = {1}\r\n", ans.n, ans.m);
            textBox_ref.Text += String.Format("метод верхней релаксации с параметром ω = {0}, применены критерии остановки по точности εмет = {1}\r\n", ans.param[0], ans.eps[0]);
            textBox_ref.Text += String.Format("и по числу итераций Nmax = {0}\r\n", ans.nmax[0]);
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("На решение схемы (СЛАУ) затрачено итераций N = {0} и достигнута точность итерационного метода ε(N) = {1}\r\n", ans.N[0], ans.acc[0]);
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("Схема (СЛАУ) решена с невязкой ||R(N)|| = {0} использована {1}\r\n", ans.R[0], "норма max");
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("Для контроля точности решения использована сетка с половинным шагом\r\n");
            textBox_ref.Text += String.Format("метод верхней релаксации с параметром ω2 = {0},\r\n", ans.param[0]);
            textBox_ref.Text += String.Format("применены критерии остановки по точности εмет-2 = {0} и по числу итераций Nmax-2 = {1}\r\n", ans.eps[1], ans.nmax[1]);
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("На решение задачи (СЛАУ) затрачено итераций N2 = {0} и достигнута точность итерационного метода ε(N2) = {1}\r\n", ans.N[1], ans.acc[1]);
            textBox_ref.Text += String.Format("Схема (СЛАУ) на сетке с половинным шагом решена с невязкой ||R(N2)|| = {0} использована норма {1}\r\n", ans.R[1], "норма max");
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("Основная задача должна быть решена с точностью не хуже чем ε = 0.5⋅10 –6;\r\n");
            textBox_ref.Text += String.Format("задача решена с точностью ε2 = {0}\r\n", ans.err);
            textBox_ref.Text += String.Format("\r\n");
            textBox_ref.Text += String.Format("Максимальное отклонение численных решений на основной сетке и сетке с\r\n");
            textBox_ref.Text += String.Format("половинным шагом наблюдается в узле x = {0}; y = {1}\r\n", ans.x_err, ans.y_err);
            textBox_ref.Text += String.Format("В качестве начального приближения на основной сетке использовано {0}, \r\n", "интерполяция по x");
            textBox_ref.Text += String.Format("на сетке с половинным шагом использовано {0}\r\n", "интерполяция по x");
            textBox_ref.Text += String.Format("\r\nВремя работы: {0} сек", ans.time);
        }

        private void ShowTable(DataGridView gridView, double start_x, double step_x, double start_y, double step_y, List<List<double>> arr, int n, int m) {
            n = Math.Min(n, 400);
            m = Math.Min(m, 400);

            gridView.ColumnCount = n + 2;
            gridView.RowCount = m + 2;

            var col = gridView.Columns[0];
            col.Name = String.Format("x[i]");
            col.ValueType = typeof(string);
            for (int i = 1; i <= n + 1; ++i) {
                col = gridView.Columns[i];
                col.Name = String.Format("x[{0}]", i);
                col.ValueType = typeof(string);
            }
            gridView.Rows[0].HeaderCell.Value = String.Format("y[j]");
            for (int j = 0; j <= m; ++j)
            {
                gridView.Rows[j + 1].HeaderCell.Value = String.Format("y[{0}]", j);
            }

            gridView[0, 0].Value = String.Format("Y\\X");

            for (int i = 0; i <= n; ++i)
                gridView[i + 1, 0].Value = String.Format("{0:f8}", start_x + i * step_x);

            for (int j = 0; j <= m; ++j)
                gridView[0, j + 1].Value = String.Format("{0:f8}", start_y + j * step_y);

            for (int j = 0; j <= m; ++j)
                for (int i = 0; i <= n; ++i)
                    gridView[i + 1, j + 1].Value = String.Format("{0:f8}", arr[j][i]);
        }

        private void ShowLastAnswer() {
            FileStream fs = File.Open(path_output, FileMode.Open, FileAccess.Read);
            StreamReader reader = new StreamReader(fs);
            if (reader == null) throw new Exception("reader is null");
            answer ans = JsonConvert.DeserializeObject<answer>(reader.ReadToEnd());
            reader.Close();
            ShowTable(dataGridView_u1, ans.a, (ans.b - ans.a) / ans.n, ans.c, (ans.d - ans.c) / ans.m, ans.arr_u[0], ans.n, ans.m);
            if (ans.test)
            {
                ShowTable(dataGridView_u2, ans.a, (ans.b - ans.a) / ans.n, ans.c, (ans.d - ans.c) / ans.m, ans.arr_u[1], ans.n, ans.m);
            }
            else
            {
                ShowTable(dataGridView_u2, ans.a, (ans.b - ans.a) / (ans.n * 2), ans.c, (ans.d - ans.c) / (ans.m * 2), ans.arr_u[1], ans.n * 2, ans.m * 2);
            }
            ShowTable(dataGridView_abs_err, ans.a, (ans.b - ans.a) / (ans.n), ans.c, (ans.d - ans.c) / (ans.m), ans.arr_err, ans.n, ans.m);

            if (ans.test)
            {
                print_ref_for_test_task(ans);
            }
            else
            {
                print_ref_for_basic_task(ans);
            }
        }

        private void backgroundWorker1_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            if (!e.Cancelled) {
                ShowLastAnswer();
            }
        }

        private void backgroundWorker1_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            progressBar1.Value = e.ProgressPercentage;
        }

        private void button2_Click(object sender, EventArgs e)
        {
            backgroundWorker1.CancelAsync();
        }

        private void Form1_FormClosing(object sender, FormClosingEventArgs e)
        {
            backgroundWorker1.CancelAsync();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            comboBox2.SelectedIndex = 0;
            comboBox1.SelectedIndex = 0;

        }

        private void comboBox2_SelectedIndexChanged(object sender, EventArgs e)
        {
            switch (comboBox2.SelectedIndex)
            {
                case 0:
                textBox_Nmax2.Enabled = false;
                textBox_eps2.Enabled = false;
                    break;
                case 1:
                textBox_Nmax2.Enabled = true;
                textBox_eps2.Enabled = true;
                    break;
            }
        }

        void ShowGraph(int var, bool test) {
            if (File.Exists(path_graph))
            {
                ProcessStartInfo info = new ProcessStartInfo();
                info.FileName = "python";
                info.WindowStyle = ProcessWindowStyle.Hidden;
                info.Arguments = String.Format(path_graph + " {0} {1}", test, var);
                Process p = Process.Start(info);
                p.WaitForExit();
            }
            else {
                MessageBox.Show(path_graph + " not found");
            }
        }

        private void button3_Click(object sender, EventArgs e)
        {
            ShowGraph(1, data.test);
        }

        private void button4_Click(object sender, EventArgs e)
        {
            ShowGraph(2, data.test);
        }

        private void button5_Click(object sender, EventArgs e)
        {
            ShowGraph(3, data.test);
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            
        }

        private void button6_Click(object sender, EventArgs e)
        {
            ShowLastAnswer();
        }
    }

    class global_data
    {
        public int numberVariant;
        public int numberMethod;
        public bool test;
        // stop condition
        public List<int> nmax;
        public List<double> eps;
        // area
        public int n, m;
        public double a, b, c, d;
        public string area;
        public int area_x, area_y;
        // method
        public List<double> param;

        public global_data() {
            nmax = new List<int>();
            eps = new List<double>();
            param = new List<double>();
        }
    }

    class answer
    {
        public double a, b, c, d;
        public int n, m;
        public double err;
        public double x_err, y_err;
        public bool test;
        public List<double> N;
        public List<double> nmax;
        public List<double> eps;
        public List<double> acc;
        public List<double> R;
        public List<double> param;
        public List<List<List<double>>> arr_u;
        public List<List<double>> arr_err;
        public double time;
    }
}
