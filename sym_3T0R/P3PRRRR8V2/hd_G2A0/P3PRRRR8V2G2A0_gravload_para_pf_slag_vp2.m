% Calculate Gravitation load for parallel robot
% P3PRRRR8V2G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:49:35
% EndTime: 2020-08-06 17:49:37
% DurationCPUTime: 1.57s
% Computational Cost: add. (927->170), mult. (1833->324), div. (36->4), fcn. (1527->22), ass. (0->125)
t933 = legFrame(3,2);
t919 = sin(t933);
t922 = cos(t933);
t904 = t919 * g(1) + t922 * g(2);
t930 = sin(pkin(4));
t932 = cos(pkin(4));
t907 = t922 * g(1) - t919 * g(2);
t929 = sin(pkin(8));
t918 = g(3) * t929;
t931 = cos(pkin(8));
t964 = t907 * t931 - t918;
t1027 = t904 * t930 + t964 * t932;
t937 = sin(qJ(2,3));
t943 = cos(qJ(2,3));
t1011 = g(3) * t931;
t965 = t907 * t929 + t1011;
t953 = t1027 * t937 + t965 * t943;
t934 = legFrame(2,2);
t920 = sin(t934);
t923 = cos(t934);
t905 = t920 * g(1) + t923 * g(2);
t908 = t923 * g(1) - t920 * g(2);
t962 = t908 * t931 - t918;
t1028 = t905 * t930 + t962 * t932;
t939 = sin(qJ(2,2));
t945 = cos(qJ(2,2));
t963 = t908 * t929 + t1011;
t952 = t1028 * t939 + t963 * t945;
t935 = legFrame(1,2);
t921 = sin(t935);
t924 = cos(t935);
t906 = t921 * g(1) + t924 * g(2);
t909 = t924 * g(1) - t921 * g(2);
t960 = t909 * t931 - t918;
t1029 = t906 * t930 + t960 * t932;
t941 = sin(qJ(2,1));
t947 = cos(qJ(2,1));
t961 = t909 * t929 + t1011;
t951 = t1029 * t941 + t961 * t947;
t936 = sin(qJ(3,3));
t1020 = pkin(2) * t936;
t938 = sin(qJ(3,2));
t1019 = pkin(2) * t938;
t940 = sin(qJ(3,1));
t1018 = pkin(2) * t940;
t942 = cos(qJ(3,3));
t1017 = pkin(3) * t942 ^ 2;
t944 = cos(qJ(3,2));
t1016 = pkin(3) * t944 ^ 2;
t946 = cos(qJ(3,1));
t1015 = pkin(3) * t946 ^ 2;
t1014 = pkin(3) * t942;
t1013 = pkin(3) * t944;
t1012 = pkin(3) * t946;
t1010 = m(3) * pkin(2) + mrSges(2,1);
t1001 = t904 * t932;
t950 = 0.1e1 / pkin(3);
t988 = t930 * t931;
t991 = mrSges(3,2) * t918 * t930;
t1009 = (((t964 * t930 - t1001) * mrSges(3,1) + t953 * mrSges(3,2)) * t942 + t936 * (t991 + (-t907 * t988 + t1001) * mrSges(3,2) + t953 * mrSges(3,1))) * t950;
t998 = t905 * t932;
t1008 = (((t962 * t930 - t998) * mrSges(3,1) + t952 * mrSges(3,2)) * t944 + t938 * (t991 + (-t908 * t988 + t998) * mrSges(3,2) + t952 * mrSges(3,1))) * t950;
t995 = t906 * t932;
t1007 = (((t960 * t930 - t995) * mrSges(3,1) + t951 * mrSges(3,2)) * t946 + t940 * (t991 + (-t909 * t988 + t995) * mrSges(3,2) + t951 * mrSges(3,1))) * t950;
t917 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t874 = -t953 * t917 + (-t1027 * t943 + t937 * t965) * (t942 * mrSges(3,1) - mrSges(3,2) * t936 + t1010);
t976 = t932 * t937;
t895 = t929 * t976 - t931 * t943;
t990 = t929 * t930;
t1006 = (t936 * t895 + t942 * t990) * t874;
t875 = -t952 * t917 + (-t1028 * t945 + t939 * t963) * (t944 * mrSges(3,1) - mrSges(3,2) * t938 + t1010);
t974 = t932 * t939;
t896 = t929 * t974 - t931 * t945;
t1005 = (t938 * t896 + t944 * t990) * t875;
t876 = -t951 * t917 + (-t1029 * t947 + t941 * t961) * (t946 * mrSges(3,1) - mrSges(3,2) * t940 + t1010);
t972 = t932 * t941;
t897 = t929 * t972 - t931 * t947;
t1004 = (t940 * t897 + t946 * t990) * t876;
t925 = m(1) + m(2) + m(3);
t1003 = t904 * t925;
t1000 = t905 * t925;
t997 = t906 * t925;
t989 = t929 * t932;
t987 = t930 * t936;
t986 = t930 * t937;
t985 = t930 * t938;
t984 = t930 * t939;
t983 = t930 * t940;
t982 = t930 * t941;
t981 = t931 * t932;
t980 = t931 * t942;
t979 = t931 * t944;
t978 = t931 * t946;
t977 = t932 * t936;
t975 = t932 * t938;
t973 = t932 * t940;
t971 = t932 * t943;
t970 = t932 * t945;
t969 = t932 * t947;
t948 = pkin(7) + pkin(6);
t910 = pkin(2) * t937 - t948 * t943;
t913 = pkin(2) * t943 + t937 * t948;
t968 = ((t929 * t971 + t931 * t937) * t1014 + t913 * t989 + t910 * t931) * t1009;
t911 = pkin(2) * t939 - t948 * t945;
t914 = pkin(2) * t945 + t939 * t948;
t967 = ((t929 * t970 + t931 * t939) * t1013 + t914 * t989 + t911 * t931) * t1008;
t912 = pkin(2) * t941 - t948 * t947;
t915 = pkin(2) * t947 + t941 * t948;
t966 = ((t929 * t969 + t931 * t941) * t1012 + t915 * t989 + t912 * t931) * t1007;
t956 = pkin(3) * t987 - t910 * t932;
t955 = pkin(3) * t985 - t911 * t932;
t954 = pkin(3) * t983 - t912 * t932;
t900 = t929 * t947 + t931 * t972;
t899 = t929 * t945 + t931 * t974;
t898 = t929 * t943 + t931 * t976;
t891 = pkin(3) * t973 + t912 * t930;
t890 = pkin(3) * t975 + t911 * t930;
t889 = pkin(3) * t977 + t910 * t930;
t885 = -t929 * t915 + t954 * t931;
t884 = -t929 * t914 + t955 * t931;
t883 = -t929 * t913 + t956 * t931;
t882 = 0.1e1 / (pkin(2) * t973 + t982 * t1015 + t891 * t946);
t881 = 0.1e1 / (pkin(2) * t975 + t984 * t1016 + t890 * t944);
t880 = 0.1e1 / (pkin(2) * t977 + t986 * t1017 + t889 * t942);
t1 = [-g(1) * m(4) + (-((t900 * t924 + t921 * t982) * t1015 + (-t885 * t924 + t891 * t921) * t946 + (t932 * t921 - t924 * t988) * t1018) * t997 - t924 * t1004 - t924 * t966) * t882 + (-((t899 * t923 + t920 * t984) * t1016 + (-t884 * t923 + t890 * t920) * t944 + (t932 * t920 - t923 * t988) * t1019) * t1000 - t923 * t1005 - t923 * t967) * t881 + (-((t898 * t922 + t919 * t986) * t1017 + (-t883 * t922 + t889 * t919) * t942 + (t932 * t919 - t922 * t988) * t1020) * t1003 - t922 * t1006 - t922 * t968) * t880; -g(2) * m(4) + (-(-(t900 * t921 - t924 * t982) * t1015 + (t885 * t921 + t924 * t891) * t946 + (t921 * t988 + t924 * t932) * t1018) * t997 + t921 * t1004 + t921 * t966) * t882 + (-(-(t899 * t920 - t923 * t984) * t1016 + (t884 * t920 + t923 * t890) * t944 + (t920 * t988 + t923 * t932) * t1019) * t1000 + t920 * t1005 + t920 * t967) * t881 + (-(-(t898 * t919 - t922 * t986) * t1017 + (t883 * t919 + t922 * t889) * t942 + (t919 * t988 + t922 * t932) * t1020) * t1003 + t919 * t1006 + t919 * t968) * t880; -g(3) * m(4) + (-(-t897 * t1015 + t915 * t978 + (pkin(2) * t983 + t954 * t946) * t929) * t997 + (-t940 * t900 - t930 * t978) * t876 + ((t929 * t941 - t931 * t969) * t1012 - t915 * t981 + t929 * t912) * t1007) * t882 + (-(-t896 * t1016 + t914 * t979 + (pkin(2) * t985 + t955 * t944) * t929) * t1000 + (-t938 * t899 - t930 * t979) * t875 + ((t929 * t939 - t931 * t970) * t1013 - t914 * t981 + t929 * t911) * t1008) * t881 + (-(-t895 * t1017 + t913 * t980 + (pkin(2) * t987 + t956 * t942) * t929) * t1003 + (-t936 * t898 - t930 * t980) * t874 + ((t929 * t937 - t931 * t971) * t1014 - t913 * t981 + t929 * t910) * t1009) * t880;];
taugX  = t1;
