% Calculate Gravitation load for parallel robot
% P3RRPRR12V1G1A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V1G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:01:37
% EndTime: 2020-08-06 19:01:37
% DurationCPUTime: 0.35s
% Computational Cost: add. (525->113), mult. (861->191), div. (30->6), fcn. (540->18), ass. (0->94)
t981 = m(1) * rSges(1,1);
t923 = (pkin(1) + rSges(3,1)) * m(3) + m(2) * rSges(2,1);
t980 = t923 * g(3);
t940 = sin(qJ(1,3));
t979 = t940 * pkin(4);
t942 = sin(qJ(1,2));
t978 = t942 * pkin(4);
t944 = sin(qJ(1,1));
t977 = t944 * pkin(4);
t939 = sin(qJ(2,3));
t945 = cos(qJ(2,3));
t953 = 0.1e1 / qJ(3,3);
t936 = legFrame(3,3);
t924 = sin(t936);
t927 = cos(t936);
t898 = -t924 * g(1) + t927 * g(2);
t901 = t927 * g(1) + t924 * g(2);
t946 = cos(qJ(1,3));
t961 = t898 * t940 + t901 * t946;
t976 = (-g(3) * t945 + t939 * t961) * t953;
t941 = sin(qJ(2,2));
t947 = cos(qJ(2,2));
t954 = 0.1e1 / qJ(3,2);
t937 = legFrame(2,3);
t925 = sin(t937);
t928 = cos(t937);
t899 = -t925 * g(1) + t928 * g(2);
t902 = t928 * g(1) + t925 * g(2);
t948 = cos(qJ(1,2));
t960 = t899 * t942 + t902 * t948;
t975 = (-g(3) * t947 + t941 * t960) * t954;
t943 = sin(qJ(2,1));
t949 = cos(qJ(2,1));
t955 = 0.1e1 / qJ(3,1);
t938 = legFrame(1,3);
t926 = sin(t938);
t929 = cos(t938);
t900 = -t926 * g(1) + t929 * g(2);
t903 = t929 * g(1) + t926 * g(2);
t950 = cos(qJ(1,1));
t959 = t900 * t944 + t903 * t950;
t974 = (-g(3) * t949 + t943 * t959) * t955;
t973 = t939 * qJ(3,3);
t972 = t941 * qJ(3,2);
t971 = t943 * qJ(3,1);
t952 = pkin(1) + pkin(2);
t970 = t952 * t945;
t969 = t952 * t947;
t968 = t952 * t949;
t951 = m(2) * rSges(2,2);
t917 = (-qJ(3,3) - rSges(3,3)) * m(3) + t951;
t967 = t953 * ((t961 * t917 - t980) * t945 + t939 * (t917 * g(3) + t961 * t923));
t918 = (-qJ(3,2) - rSges(3,3)) * m(3) + t951;
t966 = t954 * ((t960 * t918 - t980) * t947 + t941 * (t918 * g(3) + t960 * t923));
t919 = (-qJ(3,1) - rSges(3,3)) * m(3) + t951;
t965 = t955 * ((t959 * t919 - t980) * t949 + t943 * (t919 * g(3) + t959 * t923));
t964 = t945 * t967;
t963 = t947 * t966;
t962 = t949 * t965;
t958 = -t939 * t917 + t923 * t945 + t981;
t957 = -t941 * t918 + t923 * t947 + t981;
t956 = -t943 * t919 + t923 * t949 + t981;
t932 = t950 * pkin(4);
t931 = t948 * pkin(4);
t930 = t946 * pkin(4);
t916 = m(1) * rSges(1,2) - m(2) * rSges(2,3) - rSges(3,2) * m(3);
t915 = t968 + t971;
t914 = t969 + t972;
t913 = t970 + t973;
t912 = 0.1e1 / t915;
t911 = 0.1e1 / t914;
t910 = 0.1e1 / t913;
t909 = t950 * t971 - t977;
t908 = t944 * t971 + t932;
t907 = t948 * t972 - t978;
t906 = t942 * t972 + t931;
t905 = t946 * t973 - t979;
t904 = t940 * t973 + t930;
t897 = t926 * t950 + t929 * t944;
t896 = -t926 * t944 + t929 * t950;
t895 = t925 * t948 + t928 * t942;
t894 = -t925 * t942 + t928 * t948;
t893 = t924 * t946 + t927 * t940;
t892 = -t924 * t940 + t927 * t946;
t891 = t915 * t950 - t977;
t890 = t914 * t948 - t978;
t889 = t913 * t946 - t979;
t888 = t915 * t944 + t932;
t887 = t914 * t942 + t931;
t886 = t913 * t940 + t930;
t879 = (t916 * t950 + t956 * t944) * t903 + (t916 * t944 - t956 * t950) * t900;
t878 = (t916 * t948 + t957 * t942) * t902 + (t916 * t942 - t957 * t948) * t899;
t877 = (t916 * t946 + t958 * t940) * t901 + (t916 * t940 - t958 * t946) * t898;
t1 = [-m(4) * g(1) + (-t897 * t879 + (t896 * t968 - t926 * t908 + t909 * t929) * t962) * t912 + (-t895 * t878 + (t894 * t969 - t925 * t906 + t907 * t928) * t963) * t911 + (-t893 * t877 + (t892 * t970 - t924 * t904 + t905 * t927) * t964) * t910 + (-(-t926 * t888 + t929 * t891) * t974 - (-t925 * t887 + t928 * t890) * t975 - (-t924 * t886 + t927 * t889) * t976) * m(3); -m(4) * g(2) + (t896 * t879 + (t897 * t968 + t908 * t929 + t926 * t909) * t962) * t912 + (t894 * t878 + (t895 * t969 + t906 * t928 + t925 * t907) * t963) * t911 + (t892 * t877 + (t893 * t970 + t904 * t927 + t924 * t905) * t964) * t910 + (-(t888 * t929 + t926 * t891) * t974 - (t887 * t928 + t925 * t890) * t975 - (t886 * t927 + t924 * t889) * t976) * m(3); t939 * t967 + t941 * t966 + t943 * t965 - m(4) * g(3) + (-(-t949 * qJ(3,1) + t952 * t943) * t974 - (-t947 * qJ(3,2) + t952 * t941) * t975 - (-t945 * qJ(3,3) + t952 * t939) * t976) * m(3);];
taugX  = t1;
