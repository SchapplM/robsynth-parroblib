% Calculate Gravitation load for parallel robot
% P3RRPRR12V1G3A0
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
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V1G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:09:53
% EndTime: 2020-08-06 19:09:54
% DurationCPUTime: 0.80s
% Computational Cost: add. (591->122), mult. (948->215), div. (36->6), fcn. (534->18), ass. (0->97)
t922 = cos(qJ(2,3));
t929 = pkin(1) + pkin(2);
t916 = sin(qJ(2,3));
t950 = t916 * qJ(3,3);
t886 = t929 * t922 + t950;
t917 = sin(qJ(1,3));
t923 = cos(qJ(1,3));
t975 = t923 * pkin(4) + t886 * t917;
t968 = t917 * pkin(4);
t974 = t886 * t923 - t968;
t924 = cos(qJ(2,2));
t918 = sin(qJ(2,2));
t948 = t918 * qJ(3,2);
t887 = t929 * t924 + t948;
t919 = sin(qJ(1,2));
t925 = cos(qJ(1,2));
t973 = t925 * pkin(4) + t887 * t919;
t967 = t919 * pkin(4);
t972 = t887 * t925 - t967;
t926 = cos(qJ(2,1));
t920 = sin(qJ(2,1));
t946 = t920 * qJ(3,1);
t888 = t929 * t926 + t946;
t921 = sin(qJ(1,1));
t927 = cos(qJ(1,1));
t971 = t927 * pkin(4) + t888 * t921;
t966 = t921 * pkin(4);
t970 = t888 * t927 - t966;
t969 = m(1) * rSges(1,1);
t913 = legFrame(3,2);
t900 = sin(t913);
t903 = cos(t913);
t871 = t900 * g(1) + t903 * g(2);
t928 = m(2) * rSges(2,2);
t893 = (-qJ(3,3) - rSges(3,3)) * m(3) + t928;
t899 = (pkin(1) + rSges(3,1)) * m(3) + m(2) * rSges(2,1);
t930 = 0.1e1 / qJ(3,3);
t874 = t903 * g(1) - t900 * g(2);
t938 = -g(3) * t917 + t874 * t923;
t962 = ((-t899 * t871 + t938 * t893) * t922 + (t893 * t871 + t938 * t899) * t916) * t930;
t914 = legFrame(2,2);
t901 = sin(t914);
t904 = cos(t914);
t872 = t901 * g(1) + t904 * g(2);
t894 = (-qJ(3,2) - rSges(3,3)) * m(3) + t928;
t931 = 0.1e1 / qJ(3,2);
t875 = t904 * g(1) - t901 * g(2);
t937 = -g(3) * t919 + t875 * t925;
t961 = ((-t899 * t872 + t937 * t894) * t924 + (t894 * t872 + t937 * t899) * t918) * t931;
t915 = legFrame(1,2);
t902 = sin(t915);
t905 = cos(t915);
t873 = t902 * g(1) + t905 * g(2);
t895 = (-qJ(3,1) - rSges(3,3)) * m(3) + t928;
t932 = 0.1e1 / qJ(3,1);
t876 = t905 * g(1) - t902 * g(2);
t936 = -g(3) * t921 + t876 * t927;
t960 = ((-t899 * t873 + t936 * t895) * t926 + (t895 * t873 + t936 * t899) * t920) * t932;
t959 = (-t871 * t922 + t938 * t916) * t930;
t958 = (-t872 * t924 + t937 * t918) * t931;
t957 = (-t873 * t926 + t936 * t920) * t932;
t956 = t900 * qJ(3,3);
t955 = t901 * qJ(3,2);
t954 = t902 * qJ(3,1);
t953 = t903 * qJ(3,3);
t952 = t904 * qJ(3,2);
t951 = t905 * qJ(3,1);
t892 = m(1) * rSges(1,2) - m(2) * rSges(2,3) - rSges(3,2) * m(3);
t906 = g(3) * t969;
t935 = -t893 * t916 + t899 * t922;
t862 = t906 * t923 + (-t892 * t917 + t935 * t923) * g(3) + (t892 * t923 + (t935 + t969) * t917) * t874;
t949 = t917 * t862;
t934 = -t894 * t918 + t899 * t924;
t863 = t906 * t925 + (-t892 * t919 + t934 * t925) * g(3) + (t892 * t925 + (t934 + t969) * t919) * t875;
t947 = t919 * t863;
t933 = -t895 * t920 + t899 * t926;
t864 = t906 * t927 + (-t892 * t921 + t933 * t927) * g(3) + (t892 * t927 + (t933 + t969) * t921) * t876;
t945 = t921 * t864;
t944 = t923 * t929;
t943 = t925 * t929;
t942 = t927 * t929;
t941 = t929 * t916;
t940 = t929 * t918;
t939 = t929 * t920;
t909 = t926 ^ 2;
t908 = t924 ^ 2;
t907 = t922 ^ 2;
t885 = -t926 * qJ(3,1) + t939;
t884 = -t924 * qJ(3,2) + t940;
t883 = -t922 * qJ(3,3) + t941;
t882 = 0.1e1 / t888;
t881 = 0.1e1 / t887;
t880 = 0.1e1 / t886;
t879 = t927 * t946 - t966;
t878 = t925 * t948 - t967;
t877 = t923 * t950 - t968;
t1 = [-m(4) * g(1) + (-t905 * t945 + ((t905 * t942 - t954) * t909 + (t879 * t905 + t902 * t939) * t926 + t954) * t960) * t882 + (-t904 * t947 + ((t904 * t943 - t955) * t908 + (t878 * t904 + t901 * t940) * t924 + t955) * t961) * t881 + (-t903 * t949 + ((t903 * t944 - t956) * t907 + (t877 * t903 + t900 * t941) * t922 + t956) * t962) * t880 + (-(t885 * t902 + t970 * t905) * t957 - (t884 * t901 + t972 * t904) * t958 - (t883 * t900 + t974 * t903) * t959) * m(3); -m(4) * g(2) + (t902 * t945 + ((-t902 * t942 - t951) * t909 + (-t902 * t879 + t905 * t939) * t926 + t951) * t960) * t882 + (t901 * t947 + ((-t901 * t943 - t952) * t908 + (-t901 * t878 + t904 * t940) * t924 + t952) * t961) * t881 + (t900 * t949 + ((-t900 * t944 - t953) * t907 + (-t900 * t877 + t903 * t941) * t922 + t953) * t962) * t880 + (-(t885 * t905 - t970 * t902) * t957 - (t884 * t904 - t972 * t901) * t958 - (t883 * t903 - t974 * t900) * t959) * m(3); -m(4) * g(3) + (-t971 * t926 * t960 - t927 * t864) * t882 + (-t973 * t924 * t961 - t925 * t863) * t881 + (-t975 * t922 * t962 - t923 * t862) * t880 + (t971 * t957 + t973 * t958 + t975 * t959) * m(3);];
taugX  = t1;
