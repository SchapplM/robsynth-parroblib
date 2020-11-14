% Calculate Gravitation load for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
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
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:08
% EndTime: 2020-08-06 19:59:09
% DurationCPUTime: 0.56s
% Computational Cost: add. (564->124), mult. (1095->207), div. (30->9), fcn. (618->23), ass. (0->83)
t911 = sin(pkin(5));
t968 = cos(pkin(5));
t883 = -m(2) * rSges(2,1) + (-rSges(3,1) * t968 + rSges(3,2) * t911 - pkin(1)) * m(3);
t884 = m(2) * rSges(2,2) + (rSges(3,1) * t911 + t968 * rSges(3,2)) * m(3);
t925 = sin(qJ(2,1));
t931 = cos(qJ(2,1));
t976 = -t883 * t931 - t884 * t925;
t923 = sin(qJ(2,2));
t929 = cos(qJ(2,2));
t975 = -t883 * t929 - t884 * t923;
t921 = sin(qJ(2,3));
t927 = cos(qJ(2,3));
t974 = -t883 * t927 - t884 * t921;
t973 = m(1) * rSges(1,1);
t972 = m(3) * (qJ(3,3) + rSges(3,3));
t971 = m(3) * (qJ(3,2) + rSges(3,3));
t970 = m(3) * (qJ(3,1) + rSges(3,3));
t969 = pkin(2) * t911;
t918 = legFrame(3,2);
t901 = sin(t918);
t904 = cos(t918);
t888 = g(1) * t904 - g(2) * t901;
t915 = pkin(4) + qJ(3,3);
t908 = 0.1e1 / t915;
t922 = sin(qJ(1,3));
t928 = cos(qJ(1,3));
t967 = (g(3) * t922 - t888 * t928) * t908;
t919 = legFrame(2,2);
t902 = sin(t919);
t905 = cos(t919);
t889 = g(1) * t905 - g(2) * t902;
t916 = pkin(4) + qJ(3,2);
t909 = 0.1e1 / t916;
t924 = sin(qJ(1,2));
t930 = cos(qJ(1,2));
t966 = (g(3) * t924 - t889 * t930) * t909;
t920 = legFrame(1,2);
t903 = sin(t920);
t906 = cos(t920);
t890 = g(1) * t906 - g(2) * t903;
t917 = pkin(4) + qJ(3,1);
t910 = 0.1e1 / t917;
t926 = sin(qJ(1,1));
t932 = cos(qJ(1,1));
t965 = (g(3) * t926 - t890 * t932) * t910;
t885 = g(1) * t901 + g(2) * t904;
t940 = g(3) * t928 + t888 * t922;
t943 = pkin(1) * t927 + pkin(2) * cos(qJ(2,3) + pkin(5));
t958 = 0.1e1 / t943 * ((-t940 * t883 + t885 * t884) * t921 + (t885 * t883 + t940 * t884) * t927);
t886 = g(1) * t902 + g(2) * t905;
t939 = g(3) * t930 + t889 * t924;
t942 = pkin(1) * t929 + pkin(2) * cos(qJ(2,2) + pkin(5));
t957 = 0.1e1 / t942 * ((-t939 * t883 + t886 * t884) * t923 + (t886 * t883 + t939 * t884) * t929);
t887 = g(1) * t903 + g(2) * t906;
t938 = g(3) * t932 + t890 * t926;
t941 = pkin(1) * t931 + pkin(2) * cos(qJ(2,1) + pkin(5));
t956 = 0.1e1 / t941 * ((-t938 * t883 + t887 * t884) * t925 + (t887 * t883 + t938 * t884) * t931);
t895 = t968 * pkin(2) + pkin(1);
t955 = t895 * t922;
t954 = t895 * t924;
t953 = t895 * t926;
t896 = m(1) * rSges(1,2) - m(2) * rSges(2,3);
t894 = t896 * g(3);
t907 = g(3) * t973;
t952 = t908 * ((-g(3) * t972 + t894) * t928 + (t974 * g(3) + t907) * t922 + ((-t973 - t974) * t928 + (t896 - t972) * t922) * t888);
t951 = t909 * ((-g(3) * t971 + t894) * t930 + (t975 * g(3) + t907) * t924 + ((-t973 - t975) * t930 + (t896 - t971) * t924) * t889);
t950 = t910 * ((-g(3) * t970 + t894) * t932 + (t976 * g(3) + t907) * t926 + ((-t973 - t976) * t932 + (t896 - t970) * t926) * t890);
t949 = t922 * t969;
t948 = t924 * t969;
t947 = t926 * t969;
t937 = t895 * t927 - t921 * t969;
t946 = 0.1e1 / t937 * t952;
t936 = t895 * t929 - t923 * t969;
t945 = 0.1e1 / t936 * t951;
t935 = t895 * t931 - t925 * t969;
t944 = 0.1e1 / t935 * t950;
t882 = t895 * t925 + t931 * t969;
t881 = t895 * t923 + t929 * t969;
t880 = t895 * t921 + t927 * t969;
t873 = -t917 * t932 + t935 * t926;
t872 = -t930 * t916 + t936 * t924;
t871 = -t928 * t915 + t937 * t922;
t1 = [((t903 * t969 + t906 * t953) * t931 + (t895 * t903 - t906 * t947) * t925) * t944 + t903 * t956 + ((t902 * t969 + t905 * t954) * t929 + (t895 * t902 - t905 * t948) * t923) * t945 + t902 * t957 + ((t901 * t969 + t904 * t955) * t927 + (t895 * t901 - t904 * t949) * t921) * t946 + t901 * t958 - m(4) * g(1) + (-(t873 * t906 + t882 * t903) * t965 - (t872 * t905 + t881 * t902) * t966 - (t871 * t904 + t880 * t901) * t967) * m(3); ((-t903 * t953 + t906 * t969) * t931 + t925 * (t895 * t906 + t903 * t947)) * t944 + t906 * t956 + ((-t902 * t954 + t905 * t969) * t929 + t923 * (t895 * t905 + t902 * t948)) * t945 + t905 * t957 + ((-t901 * t955 + t904 * t969) * t927 + t921 * (t895 * t904 + t901 * t949)) * t946 + t904 * t958 - m(4) * g(2) + (-(-t873 * t903 + t882 * t906) * t965 - (-t872 * t902 + t881 * t905) * t966 - (-t871 * t901 + t880 * t904) * t967) * m(3); t928 * t952 + t930 * t951 + t932 * t950 - m(4) * g(3) + (-(t926 * t917 + t941 * t932) * t965 - (t924 * t916 + t942 * t930) * t966 - (t922 * t915 + t943 * t928) * t967) * m(3);];
taugX  = t1;
