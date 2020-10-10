% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x13]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:34
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR1V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:34:07
% EndTime: 2020-08-06 19:34:07
% DurationCPUTime: 0.62s
% Computational Cost: add. (438->108), mult. (738->210), div. (129->7), fcn. (789->18), ass. (0->107)
t904 = sin(qJ(1,3));
t910 = cos(qJ(1,3));
t900 = legFrame(3,2);
t884 = sin(t900);
t887 = cos(t900);
t929 = t887 * g(1) - t884 * g(2);
t866 = g(3) * t904 - t929 * t910;
t897 = pkin(3) + qJ(3,3);
t890 = 0.1e1 / t897;
t955 = t890 * t910;
t941 = t866 * t955;
t906 = sin(qJ(1,2));
t912 = cos(qJ(1,2));
t901 = legFrame(2,2);
t885 = sin(t901);
t888 = cos(t901);
t928 = t888 * g(1) - t885 * g(2);
t867 = g(3) * t906 - t928 * t912;
t898 = pkin(3) + qJ(3,2);
t891 = 0.1e1 / t898;
t953 = t891 * t912;
t940 = t867 * t953;
t908 = sin(qJ(1,1));
t914 = cos(qJ(1,1));
t902 = legFrame(1,2);
t886 = sin(t902);
t889 = cos(t902);
t927 = t889 * g(1) - t886 * g(2);
t868 = g(3) * t908 - t927 * t914;
t899 = pkin(3) + qJ(3,1);
t892 = 0.1e1 / t899;
t951 = t892 * t914;
t939 = t868 * t951;
t909 = cos(qJ(2,3));
t949 = t904 * t909;
t857 = g(3) * (pkin(1) * t949 - t910 * qJ(3,3)) - t929 * (t910 * t909 * pkin(1) + t904 * qJ(3,3));
t893 = 0.1e1 / t909;
t968 = t857 * t893;
t911 = cos(qJ(2,2));
t947 = t906 * t911;
t858 = g(3) * (pkin(1) * t947 - t912 * qJ(3,2)) - t928 * (t912 * t911 * pkin(1) + t906 * qJ(3,2));
t894 = 0.1e1 / t911;
t967 = t858 * t894;
t913 = cos(qJ(2,1));
t945 = t908 * t913;
t859 = g(3) * (pkin(1) * t945 - t914 * qJ(3,1)) - t927 * (t914 * t913 * pkin(1) + t908 * qJ(3,1));
t895 = 0.1e1 / t913;
t966 = t859 * t895;
t965 = t866 * t890;
t964 = t867 * t891;
t963 = t868 * t892;
t962 = t884 * t893;
t961 = t885 * t894;
t960 = t886 * t895;
t959 = t887 * t893;
t958 = t888 * t894;
t957 = t889 * t895;
t956 = t890 * t893;
t954 = t891 * t894;
t952 = t892 * t895;
t903 = sin(qJ(2,3));
t915 = pkin(1) + pkin(2);
t950 = t903 * t915;
t905 = sin(qJ(2,2));
t948 = t905 * t915;
t907 = sin(qJ(2,1));
t946 = t907 * t915;
t944 = t909 * t915;
t943 = t911 * t915;
t942 = t913 * t915;
t872 = -t884 * t949 + t903 * t887;
t938 = t872 * t956;
t873 = -t885 * t947 + t905 * t888;
t937 = t873 * t954;
t874 = -t886 * t945 + t907 * t889;
t936 = t874 * t952;
t875 = t903 * t884 + t887 * t949;
t935 = t875 * t956;
t876 = t905 * t885 + t888 * t947;
t934 = t876 * t954;
t877 = t907 * t886 + t889 * t945;
t933 = t877 * t952;
t932 = t866 * t903 * t956;
t931 = t867 * t905 * t954;
t930 = t868 * t907 * t952;
t926 = -t897 * t910 + t904 * t944;
t925 = -t898 * t912 + t906 * t943;
t924 = -t899 * t914 + t908 * t942;
t923 = g(3) * t910 + t929 * t904;
t922 = g(3) * t912 + t928 * t906;
t921 = g(3) * t914 + t927 * t908;
t920 = t921 * t951 + t922 * t953 + t923 * t955;
t881 = t884 * g(1) + t887 * g(2);
t851 = -t881 * t909 + t903 * t923;
t882 = t885 * g(1) + t888 * g(2);
t853 = -t882 * t911 + t905 * t922;
t883 = t886 * g(1) + t889 * g(2);
t855 = -t883 * t913 + t907 * t921;
t896 = 0.1e1 / t915;
t919 = (t851 * t962 + t853 * t961 + t855 * t960) * t896;
t918 = (t851 * t959 + t853 * t958 + t855 * t957) * t896;
t917 = t921 * t936 + t922 * t937 + t923 * t938;
t916 = t921 * t933 + t922 * t934 + t923 * t935;
t856 = t883 * t907 + t913 * t921;
t854 = t882 * t905 + t911 * t922;
t852 = t881 * t903 + t909 * t923;
t1 = [0, t866 * t935 + t867 * t934 + t868 * t933, t916, 0, 0, 0, 0, 0, t875 * t965 + t876 * t964 + t877 * t963 + t919, -t875 * t932 - t876 * t931 - t877 * t930 + (t852 * t962 + t854 * t961 + t856 * t960) * t896, -t916, (t877 * t966 - (t886 * t946 + t924 * t889) * t868) * t892 + (t876 * t967 - (t885 * t948 + t925 * t888) * t867) * t891 + (t875 * t968 - (t884 * t950 + t926 * t887) * t866) * t890 + pkin(1) * t919, -g(1); 0, t866 * t938 + t867 * t937 + t868 * t936, t917, 0, 0, 0, 0, 0, t872 * t965 + t873 * t964 + t874 * t963 + t918, -t872 * t932 - t873 * t931 - t874 * t930 + (t852 * t959 + t854 * t958 + t856 * t957) * t896, -t917, (t874 * t966 - (-t924 * t886 + t889 * t946) * t868) * t892 + (t873 * t967 - (-t925 * t885 + t888 * t948) * t867) * t891 + (t872 * t968 - (-t926 * t884 + t887 * t950) * t866) * t890 + pkin(1) * t918, -g(2); 0, t939 + t940 + t941, t920, 0, 0, 0, 0, 0, t909 * t941 + t911 * t940 + t913 * t939, -t903 * t941 - t905 * t940 - t907 * t939, -t920, (t914 * t859 - (t908 * t899 + t914 * t942) * t868) * t892 + (t912 * t858 - (t906 * t898 + t912 * t943) * t867) * t891 + (t910 * t857 - (t904 * t897 + t910 * t944) * t866) * t890, -g(3);];
tau_reg  = t1;
