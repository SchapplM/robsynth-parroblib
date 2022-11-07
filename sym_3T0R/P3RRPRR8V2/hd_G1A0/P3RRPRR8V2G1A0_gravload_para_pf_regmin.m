% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V2G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:12:12
% EndTime: 2022-11-07 13:12:13
% DurationCPUTime: 1.03s
% Computational Cost: add. (897->172), mult. (1158->287), div. (123->6), fcn. (1119->33), ass. (0->143)
t1878 = legFrame(3,3);
t1860 = sin(t1878);
t1863 = cos(t1878);
t1829 = t1860 * g(1) - t1863 * g(2);
t1832 = t1863 * g(1) + t1860 * g(2);
t1885 = sin(qJ(1,3));
t1891 = cos(qJ(1,3));
t1800 = t1885 * t1829 - t1832 * t1891;
t1879 = legFrame(2,3);
t1861 = sin(t1879);
t1864 = cos(t1879);
t1830 = t1861 * g(1) - t1864 * g(2);
t1833 = t1864 * g(1) + t1861 * g(2);
t1887 = sin(qJ(1,2));
t1893 = cos(qJ(1,2));
t1801 = t1887 * t1830 - t1833 * t1893;
t1880 = legFrame(1,3);
t1862 = sin(t1880);
t1865 = cos(t1880);
t1831 = t1862 * g(1) - t1865 * g(2);
t1834 = t1865 * g(1) + t1862 * g(2);
t1889 = sin(qJ(1,1));
t1895 = cos(qJ(1,1));
t1802 = t1889 * t1831 - t1834 * t1895;
t1894 = cos(qJ(2,1));
t1871 = t1894 * pkin(2);
t1877 = qJ(2,1) + pkin(7);
t1859 = cos(t1877);
t1934 = pkin(3) * t1859;
t1837 = 0.1e1 / (t1871 + t1934);
t1856 = sin(t1877);
t1888 = sin(qJ(2,1));
t1840 = t1888 * pkin(2) + pkin(3) * t1856;
t1883 = qJ(3,1) + pkin(5);
t1874 = -pkin(6) - t1883;
t1868 = 0.1e1 / t1874;
t1921 = t1840 * t1868;
t1903 = t1837 * t1921;
t1892 = cos(qJ(2,2));
t1870 = t1892 * pkin(2);
t1876 = qJ(2,2) + pkin(7);
t1858 = cos(t1876);
t1935 = pkin(3) * t1858;
t1836 = 0.1e1 / (t1870 + t1935);
t1855 = sin(t1876);
t1886 = sin(qJ(2,2));
t1839 = t1886 * pkin(2) + pkin(3) * t1855;
t1882 = qJ(3,2) + pkin(5);
t1873 = -pkin(6) - t1882;
t1867 = 0.1e1 / t1873;
t1922 = t1839 * t1867;
t1904 = t1836 * t1922;
t1890 = cos(qJ(2,3));
t1869 = t1890 * pkin(2);
t1875 = qJ(2,3) + pkin(7);
t1857 = cos(t1875);
t1936 = pkin(3) * t1857;
t1835 = 0.1e1 / (t1869 + t1936);
t1854 = sin(t1875);
t1884 = sin(qJ(2,3));
t1838 = t1884 * pkin(2) + pkin(3) * t1854;
t1881 = qJ(3,3) + pkin(5);
t1872 = -pkin(6) - t1881;
t1866 = 0.1e1 / t1872;
t1923 = t1838 * t1866;
t1905 = t1835 * t1923;
t1939 = t1800 * t1905 + t1801 * t1904 + t1802 * t1903;
t1828 = -t1862 * t1889 + t1865 * t1895;
t1924 = t1828 * t1868;
t1826 = -t1861 * t1887 + t1864 * t1893;
t1926 = t1826 * t1867;
t1901 = t1860 * t1885 - t1863 * t1891;
t1928 = t1901 * t1866;
t1938 = t1800 * t1928 - t1801 * t1926 - t1802 * t1924;
t1827 = t1862 * t1895 + t1865 * t1889;
t1925 = t1827 * t1868;
t1825 = t1861 * t1893 + t1864 * t1887;
t1927 = t1825 * t1867;
t1823 = t1860 * t1891 + t1863 * t1885;
t1929 = t1823 * t1866;
t1937 = t1800 * t1929 + t1801 * t1927 + t1802 * t1925;
t1933 = 0.2e1 * pkin(2) * pkin(3);
t1932 = 2 * pkin(1);
t1797 = t1823 * g(1) + t1901 * g(2);
t1931 = t1797 * t1866;
t1841 = t1887 * g(1) - t1893 * g(2);
t1842 = t1893 * g(1) + t1887 * g(2);
t1813 = t1841 * t1864 + t1842 * t1861;
t1930 = t1813 * t1867;
t1920 = t1854 * t1931;
t1919 = t1857 * t1931;
t1798 = -t1827 * g(1) + t1828 * g(2);
t1918 = t1798 * t1856 * t1868;
t1799 = -t1825 * g(1) + t1826 * g(2);
t1917 = t1799 * t1855 * t1867;
t1803 = t1829 * t1891 + t1832 * t1885;
t1916 = t1803 * t1923;
t1805 = t1831 * t1895 + t1834 * t1889;
t1915 = t1805 * t1921;
t1914 = t1858 * t1930;
t1913 = t1892 * t1930;
t1843 = g(1) * t1889 - g(2) * t1895;
t1844 = g(1) * t1895 + g(2) * t1889;
t1815 = t1843 * t1865 + t1844 * t1862;
t1912 = t1815 * t1859 * t1868;
t1911 = t1803 * t1929;
t1910 = t1803 * t1928;
t1804 = t1830 * t1893 + t1833 * t1887;
t1909 = t1804 * t1927;
t1908 = t1804 * t1926;
t1907 = t1805 * t1925;
t1906 = t1805 * t1924;
t1902 = t1797 * t1923 + g(3);
t1900 = pkin(2) ^ 2;
t1899 = pkin(3) ^ 2;
t1898 = 0.2e1 * qJ(2,1);
t1897 = 0.2e1 * qJ(2,2);
t1896 = 0.2e1 * qJ(2,3);
t1853 = t1871 + pkin(1);
t1852 = t1870 + pkin(1);
t1851 = t1869 + pkin(1);
t1850 = t1852 * t1893;
t1849 = t1851 * t1891;
t1848 = t1895 * t1853;
t1847 = t1889 * t1853;
t1846 = t1887 * t1852;
t1845 = t1885 * t1851;
t1822 = -t1887 * t1873 + t1850;
t1821 = -t1885 * t1872 + t1849;
t1820 = -t1889 * t1874 + t1848;
t1819 = t1895 * t1874 + t1847;
t1818 = t1893 * t1873 + t1846;
t1817 = t1891 * t1872 + t1845;
t1816 = -t1843 * t1862 + t1844 * t1865;
t1814 = -t1841 * t1861 + t1842 * t1864;
t1812 = (t1891 * g(1) + t1885 * g(2)) * t1863 - (t1885 * g(1) - t1891 * g(2)) * t1860;
t1796 = -g(3) * t1894 + t1816 * t1888;
t1795 = -g(3) * t1892 + t1814 * t1886;
t1794 = -g(3) * t1890 + t1812 * t1884;
t1793 = (-t1883 * t1895 + t1847) * t1834 + t1831 * (t1883 * t1889 + t1848);
t1792 = (-t1882 * t1893 + t1846) * t1833 + t1830 * (t1882 * t1887 + t1850);
t1791 = (-t1881 * t1891 + t1845) * t1832 + t1829 * (t1881 * t1885 + t1849);
t1 = [0, -t1906 - t1908 + t1910, -t1938, 0, 0, 0, 0, 0, -t1826 * t1913 + t1890 * t1910 - t1894 * t1906, -t1884 * t1910 + t1886 * t1908 + t1888 * t1906, -t1826 * t1914 - t1828 * t1912 + t1901 * t1919, -t1826 * t1917 - t1828 * t1918 - t1901 * t1920, t1938, -(t1828 * t1793 - (-t1862 * t1819 + t1820 * t1865 + t1828 * t1934) * t1805) * t1868 - (t1826 * t1792 - (-t1861 * t1818 + t1822 * t1864 + t1826 * t1935) * t1804) * t1867 - (-t1901 * t1791 - (-t1860 * t1817 + t1821 * t1863 - t1901 * t1936) * t1803) * t1866, -g(1); 0, -t1907 - t1909 - t1911, t1937, 0, 0, 0, 0, 0, -t1825 * t1913 - t1890 * t1911 - t1894 * t1907, t1884 * t1911 + t1886 * t1909 + t1888 * t1907, -t1823 * t1919 - t1825 * t1914 - t1827 * t1912, t1823 * t1920 - t1825 * t1917 - t1827 * t1918, -t1937, -(t1827 * t1793 - (t1819 * t1865 + t1820 * t1862 + t1827 * t1934) * t1805) * t1868 - (t1825 * t1792 - (t1818 * t1864 + t1822 * t1861 + t1825 * t1935) * t1804) * t1867 - (t1823 * t1791 - (t1817 * t1863 + t1821 * t1860 + t1823 * t1936) * t1803) * t1866, -g(2); 0, -t1803 * t1905 - t1804 * t1904 - t1805 * t1903, t1939, 0, 0, 0, 0, 0, (-t1894 * t1915 + t1796) * t1837 + (-t1839 * t1913 + t1795) * t1836 + (-t1890 * t1916 + t1794) * t1835, (t1816 * t1894 + (g(3) + t1915) * t1888) * t1837 + (t1814 * t1892 + (t1804 * t1922 + g(3)) * t1886) * t1836 + (t1812 * t1890 + (g(3) + t1916) * t1884) * t1835, (t1816 * t1856 + (-t1815 * t1921 - g(3)) * t1859) * t1837 + (t1814 * t1855 + (-t1813 * t1922 - g(3)) * t1858) * t1836 + (t1812 * t1854 - t1902 * t1857) * t1835, (t1816 * t1859 + (-t1798 * t1921 + g(3)) * t1856) * t1837 + (t1814 * t1858 + (-t1799 * t1922 + g(3)) * t1855) * t1836 + (t1812 * t1857 + t1902 * t1854) * t1835, -t1939, -t1793 * t1903 + (sin(t1898 + pkin(7)) * t1933 + t1899 * sin(0.2e1 * t1877) + t1900 * sin(t1898) + t1840 * t1932) * t1868 * t1837 * t1805 / 0.2e1 - t1792 * t1904 + (sin(t1897 + pkin(7)) * t1933 + t1899 * sin(0.2e1 * t1876) + t1900 * sin(t1897) + t1839 * t1932) * t1867 * t1836 * t1804 / 0.2e1 - t1791 * t1905 + (sin(t1896 + pkin(7)) * t1933 + t1899 * sin(0.2e1 * t1875) + t1900 * sin(t1896) + t1838 * t1932) * t1866 * t1835 * t1803 / 0.2e1 + (t1835 * t1794 + t1836 * t1795 + t1837 * t1796) * pkin(2), -g(3);];
tau_reg  = t1;
