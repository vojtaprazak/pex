# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\GuiMinDistClustMods.ui'
#
# Created: Tue May 16 17:33:34 2017
#      by: PyQt4 UI code generator 4.11.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1365, 580)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.tabWidget = QtGui.QTabWidget(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Carlito"))
        self.tabWidget.setFont(font)
        self.tabWidget.setToolTip(_fromUtf8(""))
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.InputOutputTab = QtGui.QWidget()
        self.InputOutputTab.setObjectName(_fromUtf8("InputOutputTab"))
        self.verticalLayout = QtGui.QVBoxLayout(self.InputOutputTab)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.verticalLayoutInput = QtGui.QVBoxLayout()
        self.verticalLayoutInput.setObjectName(_fromUtf8("verticalLayoutInput"))
        self.InputLabel = QtGui.QLabel(self.InputOutputTab)
        self.InputLabel.setObjectName(_fromUtf8("InputLabel"))
        self.verticalLayoutInput.addWidget(self.InputLabel)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.labelMOD1 = QtGui.QLabel(self.InputOutputTab)
        self.labelMOD1.setObjectName(_fromUtf8("labelMOD1"))
        self.horizontalLayout.addWidget(self.labelMOD1)
        self.lineEditMOD1 = QtGui.QLineEdit(self.InputOutputTab)
        self.lineEditMOD1.setObjectName(_fromUtf8("lineEditMOD1"))
        self.horizontalLayout.addWidget(self.lineEditMOD1)
        self.BrowseMOD1 = QtGui.QPushButton(self.InputOutputTab)
        self.BrowseMOD1.setObjectName(_fromUtf8("BrowseMOD1"))
        self.horizontalLayout.addWidget(self.BrowseMOD1)
        self.verticalLayoutInput.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.labelCSV1 = QtGui.QLabel(self.InputOutputTab)
        self.labelCSV1.setObjectName(_fromUtf8("labelCSV1"))
        self.horizontalLayout_2.addWidget(self.labelCSV1)
        self.lineEditCSV1 = QtGui.QLineEdit(self.InputOutputTab)
        self.lineEditCSV1.setObjectName(_fromUtf8("lineEditCSV1"))
        self.horizontalLayout_2.addWidget(self.lineEditCSV1)
        self.BrowseCSV1 = QtGui.QPushButton(self.InputOutputTab)
        self.BrowseCSV1.setObjectName(_fromUtf8("BrowseCSV1"))
        self.horizontalLayout_2.addWidget(self.BrowseCSV1)
        self.verticalLayoutInput.addLayout(self.horizontalLayout_2)
        self.verticalLayout.addLayout(self.verticalLayoutInput)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.InputLabel_2 = QtGui.QLabel(self.InputOutputTab)
        self.InputLabel_2.setObjectName(_fromUtf8("InputLabel_2"))
        self.verticalLayout_2.addWidget(self.InputLabel_2)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.labelMOD2 = QtGui.QLabel(self.InputOutputTab)
        self.labelMOD2.setObjectName(_fromUtf8("labelMOD2"))
        self.horizontalLayout_3.addWidget(self.labelMOD2)
        self.lineEditMOD2 = QtGui.QLineEdit(self.InputOutputTab)
        self.lineEditMOD2.setObjectName(_fromUtf8("lineEditMOD2"))
        self.horizontalLayout_3.addWidget(self.lineEditMOD2)
        self.BrowseMOD2 = QtGui.QPushButton(self.InputOutputTab)
        self.BrowseMOD2.setObjectName(_fromUtf8("BrowseMOD2"))
        self.horizontalLayout_3.addWidget(self.BrowseMOD2)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.labelCSV2 = QtGui.QLabel(self.InputOutputTab)
        self.labelCSV2.setObjectName(_fromUtf8("labelCSV2"))
        self.horizontalLayout_4.addWidget(self.labelCSV2)
        self.lineEditCSV2 = QtGui.QLineEdit(self.InputOutputTab)
        self.lineEditCSV2.setObjectName(_fromUtf8("lineEditCSV2"))
        self.horizontalLayout_4.addWidget(self.lineEditCSV2)
        self.BrowseCSV2 = QtGui.QPushButton(self.InputOutputTab)
        self.BrowseCSV2.setObjectName(_fromUtf8("BrowseCSV2"))
        self.horizontalLayout_4.addWidget(self.BrowseCSV2)
        self.verticalLayout_2.addLayout(self.horizontalLayout_4)
        self.verticalLayout.addLayout(self.verticalLayout_2)
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setObjectName(_fromUtf8("horizontalLayout_7"))
        self.labelBins = QtGui.QLabel(self.InputOutputTab)
        self.labelBins.setObjectName(_fromUtf8("labelBins"))
        self.horizontalLayout_7.addWidget(self.labelBins)
        self.spinBoxNumBins = QtGui.QSpinBox(self.InputOutputTab)
        self.spinBoxNumBins.setMinimum(2)
        self.spinBoxNumBins.setMaximum(300)
        self.spinBoxNumBins.setProperty("value", 20)
        self.spinBoxNumBins.setObjectName(_fromUtf8("spinBoxNumBins"))
        self.horizontalLayout_7.addWidget(self.spinBoxNumBins)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_7.addItem(spacerItem)
        self.labelRadius = QtGui.QLabel(self.InputOutputTab)
        self.labelRadius.setObjectName(_fromUtf8("labelRadius"))
        self.horizontalLayout_7.addWidget(self.labelRadius)
        self.lineEditRadius = QtGui.QLineEdit(self.InputOutputTab)
        self.lineEditRadius.setObjectName(_fromUtf8("lineEditRadius"))
        self.horizontalLayout_7.addWidget(self.lineEditRadius)
        self.verticalLayout.addLayout(self.horizontalLayout_7)
        self.widget_2 = QtGui.QWidget(self.InputOutputTab)
        self.widget_2.setObjectName(_fromUtf8("widget_2"))
        self.verticalLayout.addWidget(self.widget_2)
        self.StartButton = QtGui.QPushButton(self.InputOutputTab)
        self.StartButton.setObjectName(_fromUtf8("StartButton"))
        self.verticalLayout.addWidget(self.StartButton)
        spacerItem1 = QtGui.QSpacerItem(17, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        self.verticalLayout.addItem(spacerItem1)
        self.line = QtGui.QFrame(self.InputOutputTab)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.verticalLayout.addWidget(self.line)
        self.verticalLayout_5 = QtGui.QVBoxLayout()
        self.verticalLayout_5.setObjectName(_fromUtf8("verticalLayout_5"))
        self.verticalLayout_6 = QtGui.QVBoxLayout()
        self.verticalLayout_6.setObjectName(_fromUtf8("verticalLayout_6"))
        self.label_Output = QtGui.QLabel(self.InputOutputTab)
        self.label_Output.setObjectName(_fromUtf8("label_Output"))
        self.verticalLayout_6.addWidget(self.label_Output)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.label_OutMOD = QtGui.QLabel(self.InputOutputTab)
        self.label_OutMOD.setObjectName(_fromUtf8("label_OutMOD"))
        self.horizontalLayout_5.addWidget(self.label_OutMOD)
        self.lineEdit_OutMOD = QtGui.QLineEdit(self.InputOutputTab)
        self.lineEdit_OutMOD.setObjectName(_fromUtf8("lineEdit_OutMOD"))
        self.horizontalLayout_5.addWidget(self.lineEdit_OutMOD)
        self.BrowseOutMOD = QtGui.QPushButton(self.InputOutputTab)
        self.BrowseOutMOD.setObjectName(_fromUtf8("BrowseOutMOD"))
        self.horizontalLayout_5.addWidget(self.BrowseOutMOD)
        self.verticalLayout_6.addLayout(self.horizontalLayout_5)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.label_OutCSV = QtGui.QLabel(self.InputOutputTab)
        self.label_OutCSV.setObjectName(_fromUtf8("label_OutCSV"))
        self.horizontalLayout_6.addWidget(self.label_OutCSV)
        self.lineEdit_OutCSV = QtGui.QLineEdit(self.InputOutputTab)
        self.lineEdit_OutCSV.setObjectName(_fromUtf8("lineEdit_OutCSV"))
        self.horizontalLayout_6.addWidget(self.lineEdit_OutCSV)
        self.BrowseOutCSV = QtGui.QPushButton(self.InputOutputTab)
        self.BrowseOutCSV.setObjectName(_fromUtf8("BrowseOutCSV"))
        self.horizontalLayout_6.addWidget(self.BrowseOutCSV)
        self.verticalLayout_6.addLayout(self.horizontalLayout_6)
        self.verticalLayout_5.addLayout(self.verticalLayout_6)
        self.pushButton_SAVE = QtGui.QPushButton(self.InputOutputTab)
        self.pushButton_SAVE.setEnabled(False)
        self.pushButton_SAVE.setObjectName(_fromUtf8("pushButton_SAVE"))
        self.verticalLayout_5.addWidget(self.pushButton_SAVE)
        self.verticalLayout.addLayout(self.verticalLayout_5)
        spacerItem2 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem2)
        self.MytextBrowser = QtGui.QTextBrowser(self.InputOutputTab)
        self.MytextBrowser.setObjectName(_fromUtf8("MytextBrowser"))
        self.verticalLayout.addWidget(self.MytextBrowser)
        self.tabWidget.addTab(self.InputOutputTab, _fromUtf8(""))
        self.GraficTab = QtGui.QWidget()
        self.GraficTab.setObjectName(_fromUtf8("GraficTab"))
        self.verticalLayout_10 = QtGui.QVBoxLayout(self.GraficTab)
        self.verticalLayout_10.setObjectName(_fromUtf8("verticalLayout_10"))
        self.verticalLayout_9 = QtGui.QVBoxLayout()
        self.verticalLayout_9.setObjectName(_fromUtf8("verticalLayout_9"))
        self.MyGraficwidget = QtGui.QWidget(self.GraficTab)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.MyGraficwidget.sizePolicy().hasHeightForWidth())
        self.MyGraficwidget.setSizePolicy(sizePolicy)
        self.MyGraficwidget.setObjectName(_fromUtf8("MyGraficwidget"))
        self.MyGraficLayout = QtGui.QVBoxLayout(self.MyGraficwidget)
        self.MyGraficLayout.setMargin(0)
        self.MyGraficLayout.setObjectName(_fromUtf8("MyGraficLayout"))
        self.verticalLayout_9.addWidget(self.MyGraficwidget)
        self.horizontalLayout_15 = QtGui.QHBoxLayout()
        self.horizontalLayout_15.setObjectName(_fromUtf8("horizontalLayout_15"))
        self.label = QtGui.QLabel(self.GraficTab)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout_15.addWidget(self.label)
        self.doubleSpinBox_LimitLow = QtGui.QDoubleSpinBox(self.GraficTab)
        self.doubleSpinBox_LimitLow.setProperty("value", 1.0)
        self.doubleSpinBox_LimitLow.setObjectName(_fromUtf8("doubleSpinBox_LimitLow"))
        self.horizontalLayout_15.addWidget(self.doubleSpinBox_LimitLow)
        self.label_3 = QtGui.QLabel(self.GraficTab)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.horizontalLayout_15.addWidget(self.label_3)
        self.doubleSpinBox_LimitHigh = QtGui.QDoubleSpinBox(self.GraficTab)
        self.doubleSpinBox_LimitHigh.setProperty("value", 6.0)
        self.doubleSpinBox_LimitHigh.setObjectName(_fromUtf8("doubleSpinBox_LimitHigh"))
        self.horizontalLayout_15.addWidget(self.doubleSpinBox_LimitHigh)
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_15.addItem(spacerItem3)
        self.checkBox = QtGui.QCheckBox(self.GraficTab)
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.horizontalLayout_15.addWidget(self.checkBox)
        self.verticalLayout_9.addLayout(self.horizontalLayout_15)
        self.verticalLayout_10.addLayout(self.verticalLayout_9)
        self.tabWidget.addTab(self.GraficTab, _fromUtf8(""))
        self.ListTab = QtGui.QWidget()
        self.ListTab.setObjectName(_fromUtf8("ListTab"))
        self.gridLayout_4 = QtGui.QGridLayout(self.ListTab)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.tableView = QtGui.QTableView(self.ListTab)
        self.tableView.setObjectName(_fromUtf8("tableView"))
        self.gridLayout_4.addWidget(self.tableView, 0, 0, 1, 1)
        self.tableView_2 = QtGui.QTableView(self.ListTab)
        self.tableView_2.setObjectName(_fromUtf8("tableView_2"))
        self.gridLayout_4.addWidget(self.tableView_2, 1, 0, 1, 1)
        self.tabWidget.addTab(self.ListTab, _fromUtf8(""))
        self.tab_cluster_input = QtGui.QWidget()
        self.tab_cluster_input.setObjectName(_fromUtf8("tab_cluster_input"))
        self.gridLayout_2 = QtGui.QGridLayout(self.tab_cluster_input)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.horizontalLayout_17 = QtGui.QHBoxLayout()
        self.horizontalLayout_17.setObjectName(_fromUtf8("horizontalLayout_17"))
        self.label_2 = QtGui.QLabel(self.tab_cluster_input)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.horizontalLayout_17.addWidget(self.label_2)
        self.doubleSpinBox_min_limit_clust = QtGui.QDoubleSpinBox(self.tab_cluster_input)
        self.doubleSpinBox_min_limit_clust.setProperty("value", 7.5)
        self.doubleSpinBox_min_limit_clust.setObjectName(_fromUtf8("doubleSpinBox_min_limit_clust"))
        self.horizontalLayout_17.addWidget(self.doubleSpinBox_min_limit_clust)
        spacerItem4 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_17.addItem(spacerItem4)
        self.label_7 = QtGui.QLabel(self.tab_cluster_input)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.horizontalLayout_17.addWidget(self.label_7)
        self.doubleSpinBox_max_limit_cluster = QtGui.QDoubleSpinBox(self.tab_cluster_input)
        self.doubleSpinBox_max_limit_cluster.setProperty("value", 9.5)
        self.doubleSpinBox_max_limit_cluster.setObjectName(_fromUtf8("doubleSpinBox_max_limit_cluster"))
        self.horizontalLayout_17.addWidget(self.doubleSpinBox_max_limit_cluster)
        spacerItem5 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_17.addItem(spacerItem5)
        self.checkBox_Cluster_nearest_only = QtGui.QCheckBox(self.tab_cluster_input)
        self.checkBox_Cluster_nearest_only.setObjectName(_fromUtf8("checkBox_Cluster_nearest_only"))
        self.horizontalLayout_17.addWidget(self.checkBox_Cluster_nearest_only)
        spacerItem6 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_17.addItem(spacerItem6)
        self.verticalLayout_3.addLayout(self.horizontalLayout_17)
        self.horizontalLayout_8 = QtGui.QHBoxLayout()
        self.horizontalLayout_8.setObjectName(_fromUtf8("horizontalLayout_8"))
        self.label_4 = QtGui.QLabel(self.tab_cluster_input)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.horizontalLayout_8.addWidget(self.label_4)
        self.comboBox_Cluster = QtGui.QComboBox(self.tab_cluster_input)
        self.comboBox_Cluster.setObjectName(_fromUtf8("comboBox_Cluster"))
        self.comboBox_Cluster.addItem(_fromUtf8(""))
        self.comboBox_Cluster.addItem(_fromUtf8(""))
        self.comboBox_Cluster.addItem(_fromUtf8(""))
        self.comboBox_Cluster.addItem(_fromUtf8(""))
        self.comboBox_Cluster.addItem(_fromUtf8(""))
        self.horizontalLayout_8.addWidget(self.comboBox_Cluster)
        self.verticalLayout_3.addLayout(self.horizontalLayout_8)
        self.horizontalLayout_9 = QtGui.QHBoxLayout()
        self.horizontalLayout_9.setContentsMargins(-1, -1, -1, 0)
        self.horizontalLayout_9.setObjectName(_fromUtf8("horizontalLayout_9"))
        self.label_5 = QtGui.QLabel(self.tab_cluster_input)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.horizontalLayout_9.addWidget(self.label_5)
        self.comboBox_cluster_what = QtGui.QComboBox(self.tab_cluster_input)
        self.comboBox_cluster_what.setObjectName(_fromUtf8("comboBox_cluster_what"))
        self.comboBox_cluster_what.addItem(_fromUtf8(""))
        self.comboBox_cluster_what.addItem(_fromUtf8(""))
        self.horizontalLayout_9.addWidget(self.comboBox_cluster_what)
        self.verticalLayout_3.addLayout(self.horizontalLayout_9)
        self.checkBox_Cluster_inv_rot = QtGui.QCheckBox(self.tab_cluster_input)
        self.checkBox_Cluster_inv_rot.setObjectName(_fromUtf8("checkBox_Cluster_inv_rot"))
        self.verticalLayout_3.addWidget(self.checkBox_Cluster_inv_rot)
        self.pushButton_StartCluster = QtGui.QPushButton(self.tab_cluster_input)
        self.pushButton_StartCluster.setObjectName(_fromUtf8("pushButton_StartCluster"))
        self.verticalLayout_3.addWidget(self.pushButton_StartCluster)
        self.horizontalLayout_19 = QtGui.QHBoxLayout()
        self.horizontalLayout_19.setObjectName(_fromUtf8("horizontalLayout_19"))
        self.radioButton_cluster_number = QtGui.QRadioButton(self.tab_cluster_input)
        self.radioButton_cluster_number.setChecked(True)
        self.radioButton_cluster_number.setObjectName(_fromUtf8("radioButton_cluster_number"))
        self.horizontalLayout_19.addWidget(self.radioButton_cluster_number)
        self.spinBox_number_cluster = QtGui.QSpinBox(self.tab_cluster_input)
        self.spinBox_number_cluster.setProperty("value", 4)
        self.spinBox_number_cluster.setObjectName(_fromUtf8("spinBox_number_cluster"))
        self.horizontalLayout_19.addWidget(self.spinBox_number_cluster)
        spacerItem7 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_19.addItem(spacerItem7)
        self.radioButton_cluster_distance = QtGui.QRadioButton(self.tab_cluster_input)
        self.radioButton_cluster_distance.setObjectName(_fromUtf8("radioButton_cluster_distance"))
        self.horizontalLayout_19.addWidget(self.radioButton_cluster_distance)
        self.doubleSpinBox_Cluster_MaxDistance = QtGui.QDoubleSpinBox(self.tab_cluster_input)
        self.doubleSpinBox_Cluster_MaxDistance.setSingleStep(0.01)
        self.doubleSpinBox_Cluster_MaxDistance.setProperty("value", 0.1)
        self.doubleSpinBox_Cluster_MaxDistance.setObjectName(_fromUtf8("doubleSpinBox_Cluster_MaxDistance"))
        self.horizontalLayout_19.addWidget(self.doubleSpinBox_Cluster_MaxDistance)
        spacerItem8 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_19.addItem(spacerItem8)
        self.verticalLayout_3.addLayout(self.horizontalLayout_19)
        self.horizontalLayout_10 = QtGui.QHBoxLayout()
        self.horizontalLayout_10.setContentsMargins(-1, -1, -1, 0)
        self.horizontalLayout_10.setObjectName(_fromUtf8("horizontalLayout_10"))
        self.label_6 = QtGui.QLabel(self.tab_cluster_input)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.horizontalLayout_10.addWidget(self.label_6)
        self.lineEdit_ClusterOut = QtGui.QLineEdit(self.tab_cluster_input)
        self.lineEdit_ClusterOut.setObjectName(_fromUtf8("lineEdit_ClusterOut"))
        self.horizontalLayout_10.addWidget(self.lineEdit_ClusterOut)
        self.pushButton_Cluster_browse = QtGui.QPushButton(self.tab_cluster_input)
        self.pushButton_Cluster_browse.setObjectName(_fromUtf8("pushButton_Cluster_browse"))
        self.horizontalLayout_10.addWidget(self.pushButton_Cluster_browse)
        self.verticalLayout_3.addLayout(self.horizontalLayout_10)
        self.pushButton_SaveCluster = QtGui.QPushButton(self.tab_cluster_input)
        self.pushButton_SaveCluster.setObjectName(_fromUtf8("pushButton_SaveCluster"))
        self.verticalLayout_3.addWidget(self.pushButton_SaveCluster)
        self.textBrowser_Cluster = QtGui.QTextBrowser(self.tab_cluster_input)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.textBrowser_Cluster.sizePolicy().hasHeightForWidth())
        self.textBrowser_Cluster.setSizePolicy(sizePolicy)
        self.textBrowser_Cluster.setObjectName(_fromUtf8("textBrowser_Cluster"))
        self.verticalLayout_3.addWidget(self.textBrowser_Cluster)
        self.gridLayout_2.addLayout(self.verticalLayout_3, 0, 0, 1, 1)
        self.widget_grafic_cluster = QtGui.QWidget(self.tab_cluster_input)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widget_grafic_cluster.sizePolicy().hasHeightForWidth())
        self.widget_grafic_cluster.setSizePolicy(sizePolicy)
        self.widget_grafic_cluster.setObjectName(_fromUtf8("widget_grafic_cluster"))
        self.verticalLayout_ClusterGrafic = QtGui.QVBoxLayout(self.widget_grafic_cluster)
        self.verticalLayout_ClusterGrafic.setMargin(0)
        self.verticalLayout_ClusterGrafic.setObjectName(_fromUtf8("verticalLayout_ClusterGrafic"))
        self.gridLayout_2.addWidget(self.widget_grafic_cluster, 0, 1, 1, 1)
        self.tabWidget.addTab(self.tab_cluster_input, _fromUtf8(""))
        self.gridLayout.addWidget(self.tabWidget, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionClose = QtGui.QAction(MainWindow)
        self.actionClose.setObjectName(_fromUtf8("actionClose"))

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(3)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "Minimal distances between two MOD models", None))
        self.InputLabel.setText(_translate("MainWindow", "Input Object 1", None))
        self.labelMOD1.setText(_translate("MainWindow", "MOD file", None))
        self.lineEditMOD1.setText(_translate("MainWindow", "D:/Data/ExperimentalData/Oxford/TransportMichael/G53a_b4_ves1_fromIter5_bin0.5_fromIter0_bin0.5_MOTL_Tom1_Iter2_remdup_2.0.mod", None))
        self.BrowseMOD1.setText(_translate("MainWindow", "browse", None))
        self.labelCSV1.setText(_translate("MainWindow", "CSV file", None))
        self.lineEditCSV1.setText(_translate("MainWindow", "D:/Data/ExperimentalData/Oxford/TransportMichael/G53a_b4_ves1_fromIter5_bin0.5_fromIter0_bin0.5_fromIter2_remdup2.0_MOTL_Tom1_Iter6.csv", None))
        self.BrowseCSV1.setText(_translate("MainWindow", "browse", None))
        self.InputLabel_2.setText(_translate("MainWindow", "Input Object 2", None))
        self.labelMOD2.setText(_translate("MainWindow", "MOD file", None))
        self.lineEditMOD2.setText(_translate("MainWindow", "D:/Data/ExperimentalData/Oxford/TransportMichael/run1_b4_NEC_C-capsids_fromIter7_hexons_fromIter0_remedge2.0_MOTL_Tom9_Iter6_remdup_0.0.mod", None))
        self.BrowseMOD2.setText(_translate("MainWindow", "browse", None))
        self.labelCSV2.setText(_translate("MainWindow", "CSV file", None))
        self.lineEditCSV2.setText(_translate("MainWindow", "D:/Data/ExperimentalData/Oxford/TransportMichael/run1_b4_NEC_C-capsids_fromIter7_hexons_fromIter0_remedge2.0_fromIter6_remdup0.0_MOTL_Tom9_Iter5.csv", None))
        self.BrowseCSV2.setText(_translate("MainWindow", "browse", None))
        self.labelBins.setText(_translate("MainWindow", "Number of bins", None))
        self.labelRadius.setText(_translate("MainWindow", "Neighbourhood radius", None))
        self.lineEditRadius.setText(_translate("MainWindow", "9", None))
        self.StartButton.setText(_translate("MainWindow", "Start", None))
        self.label_Output.setText(_translate("MainWindow", "Output", None))
        self.label_OutMOD.setText(_translate("MainWindow", "MOD file", None))
        self.lineEdit_OutMOD.setText(_translate("MainWindow", "dumbmemod", None))
        self.BrowseOutMOD.setText(_translate("MainWindow", "browse", None))
        self.label_OutCSV.setText(_translate("MainWindow", "CSV file", None))
        self.lineEdit_OutCSV.setText(_translate("MainWindow", "dumbmecsv", None))
        self.BrowseOutCSV.setText(_translate("MainWindow", "browse", None))
        self.pushButton_SAVE.setText(_translate("MainWindow", "Save", None))
        self.MytextBrowser.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Carlito\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:\'Cantarell\'; font-size:11pt;\">This program compares the positions of points in two model files. It computes the nearest neighbours of the first file and prints out some histograms.</span></p></body></html>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.InputOutputTab), _translate("MainWindow", "Input / Output", None))
        self.tabWidget.setTabToolTip(self.tabWidget.indexOf(self.InputOutputTab), _translate("MainWindow", "File Input and Output. You have to start here!", None))
        self.label.setText(_translate("MainWindow", "lower limit", None))
        self.label_3.setText(_translate("MainWindow", "upper limit", None))
        self.checkBox.setText(_translate("MainWindow", "show limits", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.GraficTab), _translate("MainWindow", "Grafics", None))
        self.tabWidget.setTabToolTip(self.tabWidget.indexOf(self.GraficTab), _translate("MainWindow", "Histogram of the relative distances of both mod files. ", None))
        self.tableView.setToolTip(_translate("MainWindow", "Nothing to see here (yet)! Move along. =>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.ListTab), _translate("MainWindow", "Listoutput", None))
        self.tabWidget.setTabToolTip(self.tabWidget.indexOf(self.ListTab), _translate("MainWindow", "Not yet implemented ", None))
        self.label_2.setText(_translate("MainWindow", "Lower limit", None))
        self.label_7.setText(_translate("MainWindow", "Upper Limit", None))
        self.checkBox_Cluster_nearest_only.setText(_translate("MainWindow", "Use only nearest neighbours", None))
        self.label_4.setText(_translate("MainWindow", "Method (see below for details)", None))
        self.comboBox_Cluster.setItemText(0, _translate("MainWindow", "single", None))
        self.comboBox_Cluster.setItemText(1, _translate("MainWindow", "complete", None))
        self.comboBox_Cluster.setItemText(2, _translate("MainWindow", "average", None))
        self.comboBox_Cluster.setItemText(3, _translate("MainWindow", "centroid", None))
        self.comboBox_Cluster.setItemText(4, _translate("MainWindow", "median", None))
        self.label_5.setText(_translate("MainWindow", "Cluster relative positions or orientations", None))
        self.comboBox_cluster_what.setItemText(0, _translate("MainWindow", "Orientations", None))
        self.comboBox_cluster_what.setItemText(1, _translate("MainWindow", "Positions", None))
        self.checkBox_Cluster_inv_rot.setText(_translate("MainWindow", "InvertRotation for Positions and Matrix multiplication from the right for Orientations", None))
        self.pushButton_StartCluster.setToolTip(_translate("MainWindow", "Starts the cluster algorithm. For large data sets the computation may take some time!!!", None))
        self.pushButton_StartCluster.setText(_translate("MainWindow", "Start", None))
        self.radioButton_cluster_number.setText(_translate("MainWindow", "Number of cluster", None))
        self.radioButton_cluster_distance.setText(_translate("MainWindow", "Max distance", None))
        self.label_6.setText(_translate("MainWindow", "Output file name", None))
        self.lineEdit_ClusterOut.setText(_translate("MainWindow", "mydumyfile", None))
        self.pushButton_Cluster_browse.setText(_translate("MainWindow", "browse", None))
        self.pushButton_SaveCluster.setToolTip(_translate("MainWindow", "Saves the resulting model and orentations to a csv and mod file. Names are defined in the first tab!", None))
        self.pushButton_SaveCluster.setText(_translate("MainWindow", "Save to Mod and CSV", None))
        self.textBrowser_Cluster.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Carlito\'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:10pt; font-weight:600;\">Clustering according to the orientation of both particles</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">------------------------------------------------------------------------------------------------------</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">method=’single’ assign d(u,v) = min(dist(u[i],v[j]))</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">method=’complete’ assigns d(u, v) = max(dist(u[i],v[j]))</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">method=’average’ assigns d(u,v) = sum_{ij} {d(u[i], v[j])}/ {(|u|*|v|)}</span></p>\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">method=’centroid’ assigns  dist(s,t) = ||c_s-c_t||</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">method=’median’ assigns math:</span><span style=\" font-size:8pt; font-style:italic;\">d(s,t)</span><span style=\" font-size:8pt;\"> like the </span><span style=\" font-family:\'Courier New,courier\'; font-size:8pt;\">centroid</span><span style=\" font-size:8pt;\"> method.</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">see also</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">or</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">------------------------------------------------------------------------------------------------------</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:8pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">You can choose to cluster in regard to the relative position of the second model to the oriented first model particles (nearest neighbour)</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">Alternatively only the relative orientations can also be clustered. </span></p></body></html>", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_cluster_input), _translate("MainWindow", "Clustering Parameters", None))
        self.tabWidget.setTabToolTip(self.tabWidget.indexOf(self.tab_cluster_input), _translate("MainWindow", "Here you can cluster the data. ", None))
        self.actionClose.setText(_translate("MainWindow", "Close", None))

