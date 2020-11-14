% Calculate Gravitation load for parallel robot
% P3PRRRR8V2G3A0
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
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V2G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:04:41
% EndTime: 2020-08-06 18:04:43
% DurationCPUTime: 1.41s
% Computational Cost: add. (927->176), mult. (1833->333), div. (36->4), fcn. (1527->22), ass. (0->131)
t929 = legFrame(3,2);
t915 = sin(t929);
t918 = cos(t929);
t899 = t915 * g(1) + t918 * g(2);
t926 = sin(pkin(4));
t928 = cos(pkin(4));
t902 = t918 * g(1) - t915 * g(2);
t927 = cos(pkin(8));
t914 = g(3) * t927;
t925 = sin(pkin(8));
t961 = -t902 * t925 - t914;
t1022 = t899 * t926 + t961 * t928;
t933 = sin(qJ(2,3));
t1027 = t1022 * t933;
t930 = legFrame(2,2);
t916 = sin(t930);
t919 = cos(t930);
t900 = t916 * g(1) + t919 * g(2);
t903 = t919 * g(1) - t916 * g(2);
t959 = -t903 * t925 - t914;
t1023 = t900 * t926 + t959 * t928;
t935 = sin(qJ(2,2));
t1026 = t1023 * t935;
t931 = legFrame(1,2);
t917 = sin(t931);
t920 = cos(t931);
t901 = t917 * g(1) + t920 * g(2);
t904 = t920 * g(1) - t917 * g(2);
t957 = -t904 * t925 - t914;
t1024 = t901 * t926 + t957 * t928;
t937 = sin(qJ(2,1));
t1025 = t1024 * t937;
t932 = sin(qJ(3,3));
t1018 = pkin(2) * t932;
t934 = sin(qJ(3,2));
t1017 = pkin(2) * t934;
t936 = sin(qJ(3,1));
t1016 = pkin(2) * t936;
t938 = cos(qJ(3,3));
t1015 = pkin(3) * t938 ^ 2;
t940 = cos(qJ(3,2));
t1014 = pkin(3) * t940 ^ 2;
t942 = cos(qJ(3,1));
t1013 = pkin(3) * t942 ^ 2;
t1012 = pkin(3) * t938;
t1011 = pkin(3) * t940;
t1010 = pkin(3) * t942;
t1009 = g(3) * t925;
t1008 = m(3) * pkin(2) + mrSges(2,1);
t946 = 0.1e1 / pkin(3);
t939 = cos(qJ(2,3));
t960 = t902 * t927 - t1009;
t949 = t960 * t939 + t1027;
t988 = t925 * t926;
t989 = mrSges(3,2) * t926 * t914;
t999 = t899 * t928;
t1007 = (((t961 * t926 - t999) * mrSges(3,1) + t949 * mrSges(3,2)) * t938 + t932 * (t989 + (t902 * t988 + t999) * mrSges(3,2) + t949 * mrSges(3,1))) * t946;
t941 = cos(qJ(2,2));
t958 = t903 * t927 - t1009;
t948 = t958 * t941 + t1026;
t996 = t900 * t928;
t1006 = (((t959 * t926 - t996) * mrSges(3,1) + t948 * mrSges(3,2)) * t940 + t934 * (t989 + (t903 * t988 + t996) * mrSges(3,2) + t948 * mrSges(3,1))) * t946;
t943 = cos(qJ(2,1));
t956 = t904 * t927 - t1009;
t947 = t956 * t943 + t1025;
t993 = t901 * t928;
t1005 = (((t957 * t926 - t993) * mrSges(3,1) + t947 * mrSges(3,2)) * t942 + t936 * (t989 + (t904 * t988 + t993) * mrSges(3,2) + t947 * mrSges(3,1))) * t946;
t913 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t911 = t913 * t1009;
t976 = t927 * t939;
t869 = t911 * t939 + (-t902 * t976 - t1027) * t913 + (-t1022 * t939 + t933 * t960) * (t938 * mrSges(3,1) - mrSges(3,2) * t932 + t1008);
t972 = t928 * t933;
t893 = t925 * t939 + t927 * t972;
t980 = t926 * t938;
t1004 = (t932 * t893 + t927 * t980) * t869;
t975 = t927 * t941;
t870 = t911 * t941 + (-t903 * t975 - t1026) * t913 + (-t1023 * t941 + t935 * t958) * (t940 * mrSges(3,1) - mrSges(3,2) * t934 + t1008);
t970 = t928 * t935;
t894 = t925 * t941 + t927 * t970;
t979 = t926 * t940;
t1003 = (t934 * t894 + t927 * t979) * t870;
t974 = t927 * t943;
t871 = t911 * t943 + (-t904 * t974 - t1025) * t913 + (-t1024 * t943 + t937 * t956) * (t942 * mrSges(3,1) - mrSges(3,2) * t936 + t1008);
t968 = t928 * t937;
t895 = t925 * t943 + t927 * t968;
t978 = t926 * t942;
t1002 = (t936 * t895 + t927 * t978) * t871;
t921 = m(1) + m(2) + m(3);
t1001 = t899 * t921;
t998 = t900 * t921;
t995 = t901 * t921;
t987 = t925 * t928;
t986 = t926 * t932;
t985 = t926 * t933;
t984 = t926 * t934;
t983 = t926 * t935;
t982 = t926 * t936;
t981 = t926 * t937;
t977 = t927 * t928;
t973 = t928 * t932;
t971 = t928 * t934;
t969 = t928 * t936;
t967 = t928 * t939;
t966 = t928 * t941;
t965 = t928 * t943;
t944 = pkin(7) + pkin(6);
t905 = pkin(2) * t933 - t944 * t939;
t908 = pkin(2) * t939 + t933 * t944;
t964 = ((t925 * t933 - t927 * t967) * t1012 - t908 * t977 + t905 * t925) * t1007;
t906 = pkin(2) * t935 - t944 * t941;
t909 = pkin(2) * t941 + t935 * t944;
t963 = ((t925 * t935 - t927 * t966) * t1011 - t909 * t977 + t906 * t925) * t1006;
t907 = pkin(2) * t937 - t944 * t943;
t910 = pkin(2) * t943 + t937 * t944;
t962 = ((t925 * t937 - t927 * t965) * t1010 - t910 * t977 + t907 * t925) * t1005;
t952 = pkin(3) * t986 - t905 * t928;
t951 = pkin(3) * t984 - t906 * t928;
t950 = pkin(3) * t982 - t907 * t928;
t892 = t925 * t968 - t974;
t891 = t925 * t970 - t975;
t890 = t925 * t972 - t976;
t886 = pkin(3) * t969 + t926 * t907;
t885 = pkin(3) * t971 + t926 * t906;
t884 = pkin(3) * t973 + t926 * t905;
t880 = t927 * t910 + t950 * t925;
t879 = t927 * t909 + t951 * t925;
t878 = t927 * t908 + t952 * t925;
t877 = 0.1e1 / (pkin(2) * t969 + t981 * t1013 + t886 * t942);
t876 = 0.1e1 / (pkin(2) * t971 + t983 * t1014 + t885 * t940);
t875 = 0.1e1 / (pkin(2) * t973 + t985 * t1015 + t884 * t938);
t1 = [-g(1) * m(4) + (-(-(t892 * t920 - t917 * t981) * t1013 + (t880 * t920 + t917 * t886) * t942 + (t928 * t917 + t920 * t988) * t1016) * t995 - t920 * t1002 + t920 * t962) * t877 + (-(-(t891 * t919 - t916 * t983) * t1014 + (t879 * t919 + t916 * t885) * t940 + (t928 * t916 + t919 * t988) * t1017) * t998 - t919 * t1003 + t919 * t963) * t876 + (-(-(t890 * t918 - t915 * t985) * t1015 + (t878 * t918 + t915 * t884) * t938 + (t928 * t915 + t918 * t988) * t1018) * t1001 - t918 * t1004 + t918 * t964) * t875; -g(2) * m(4) + (-((t892 * t917 + t920 * t981) * t1013 + (-t880 * t917 + t920 * t886) * t942 + (-t917 * t988 + t920 * t928) * t1016) * t995 + t917 * t1002 - t917 * t962) * t877 + (-((t891 * t916 + t919 * t983) * t1014 + (-t879 * t916 + t919 * t885) * t940 + (-t916 * t988 + t919 * t928) * t1017) * t998 + t916 * t1003 - t916 * t963) * t876 + (-((t890 * t915 + t918 * t985) * t1015 + (-t878 * t915 + t918 * t884) * t938 + (-t915 * t988 + t918 * t928) * t1018) * t1001 + t915 * t1004 - t915 * t964) * t875; -g(3) * m(4) + (-(-t895 * t1013 - t910 * t925 * t942 + (pkin(2) * t982 + t950 * t942) * t927) * t995 + (t936 * t892 + t925 * t978) * t871 + ((t925 * t965 + t927 * t937) * t1010 + t910 * t987 + t907 * t927) * t1005) * t877 + (-(-t894 * t1014 - t909 * t925 * t940 + (pkin(2) * t984 + t951 * t940) * t927) * t998 + (t934 * t891 + t925 * t979) * t870 + ((t925 * t966 + t927 * t935) * t1011 + t909 * t987 + t906 * t927) * t1006) * t876 + (-(-t893 * t1015 - t908 * t925 * t938 + (pkin(2) * t986 + t952 * t938) * t927) * t1001 + (t932 * t890 + t925 * t980) * t869 + ((t925 * t967 + t927 * t933) * t1012 + t908 * t987 + t905 * t927) * t1007) * t875;];
taugX  = t1;
