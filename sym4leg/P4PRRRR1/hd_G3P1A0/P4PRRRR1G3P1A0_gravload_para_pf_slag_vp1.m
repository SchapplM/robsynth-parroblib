% Calculate Gravitation load for parallel robot
% P4PRRRR1G3P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR1G3P1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:00:59
% EndTime: 2020-03-02 19:01:01
% DurationCPUTime: 1.35s
% Computational Cost: add. (467->120), mult. (966->233), div. (92->13), fcn. (698->26), ass. (0->102)
t963 = sin(qJ(3,1));
t969 = cos(qJ(3,1));
t1012 = -rSges(3,1) * t969 + rSges(3,2) * t963;
t961 = sin(qJ(3,2));
t967 = cos(qJ(3,2));
t1011 = -rSges(3,1) * t967 + rSges(3,2) * t961;
t959 = sin(qJ(3,3));
t965 = cos(qJ(3,3));
t1010 = -rSges(3,1) * t965 + rSges(3,2) * t959;
t951 = sin(qJ(3,4));
t953 = cos(qJ(3,4));
t1009 = -rSges(3,1) * t953 + rSges(3,2) * t951;
t1008 = rSges(3,1) * g(3);
t955 = legFrame(4,2);
t932 = sin(t955);
t936 = cos(t955);
t924 = t932 * g(1) + t936 * g(2);
t952 = sin(qJ(2,4));
t942 = 0.1e1 / t952;
t999 = t924 * t942;
t956 = legFrame(3,2);
t933 = sin(t956);
t937 = cos(t956);
t925 = t933 * g(1) + t937 * g(2);
t960 = sin(qJ(2,3));
t945 = 0.1e1 / t960;
t998 = t925 * t945;
t957 = legFrame(2,2);
t934 = sin(t957);
t938 = cos(t957);
t926 = t934 * g(1) + t938 * g(2);
t962 = sin(qJ(2,2));
t946 = 0.1e1 / t962;
t997 = t926 * t946;
t958 = legFrame(1,2);
t935 = sin(t958);
t939 = cos(t958);
t927 = t935 * g(1) + t939 * g(2);
t964 = sin(qJ(2,1));
t947 = 0.1e1 / t964;
t996 = t927 * t947;
t995 = t942 * t951;
t994 = t945 * t959;
t993 = t946 * t961;
t992 = t947 * t963;
t928 = t936 * g(1) - t932 * g(2);
t954 = cos(qJ(2,4));
t904 = ((-rSges(2,1) * t924 + rSges(2,2) * t928) * t954 + t952 * (rSges(2,1) * t928 + rSges(2,2) * t924)) * m(2) + ((-t928 * rSges(3,3) + t1009 * t924) * t954 + t952 * (-t924 * rSges(3,3) - t1009 * t928)) * m(3);
t943 = 0.1e1 / t953;
t991 = t904 * t942 * t943;
t929 = t937 * g(1) - t933 * g(2);
t966 = cos(qJ(2,3));
t905 = ((-rSges(2,1) * t925 + rSges(2,2) * t929) * t966 + t960 * (rSges(2,1) * t929 + rSges(2,2) * t925)) * m(2) + ((-t929 * rSges(3,3) + t1010 * t925) * t966 + t960 * (-t925 * rSges(3,3) - t1010 * t929)) * m(3);
t948 = 0.1e1 / t965;
t990 = t905 * t945 * t948;
t930 = t938 * g(1) - t934 * g(2);
t968 = cos(qJ(2,2));
t906 = ((-rSges(2,1) * t926 + rSges(2,2) * t930) * t968 + t962 * (rSges(2,1) * t930 + rSges(2,2) * t926)) * m(2) + ((-t930 * rSges(3,3) + t1011 * t926) * t968 + t962 * (-t926 * rSges(3,3) - t1011 * t930)) * m(3);
t949 = 0.1e1 / t967;
t989 = t906 * t946 * t949;
t931 = t939 * g(1) - t935 * g(2);
t970 = cos(qJ(2,1));
t907 = ((-rSges(2,1) * t927 + rSges(2,2) * t931) * t970 + t964 * (rSges(2,1) * t931 + rSges(2,2) * t927)) * m(2) + ((-t931 * rSges(3,3) + t1012 * t927) * t970 + t964 * (-t927 * rSges(3,3) - t1012 * t931)) * m(3);
t950 = 0.1e1 / t969;
t988 = t907 * t947 * t950;
t987 = t952 * t924 + t928 * t954;
t986 = t960 * t925 + t929 * t966;
t985 = t962 * t926 + t930 * t968;
t984 = t964 * t927 + t931 * t970;
t983 = 0.1e1 / pkin(2);
t982 = koppelP(1,1);
t981 = koppelP(2,1);
t980 = koppelP(3,1);
t979 = koppelP(4,1);
t978 = koppelP(1,2);
t977 = koppelP(2,2);
t976 = koppelP(3,2);
t975 = koppelP(4,2);
t974 = rSges(4,1);
t973 = rSges(4,2);
t972 = xP(4);
t971 = rSges(3,2) * g(3);
t944 = m(1) + m(2) + m(3);
t941 = cos(t972);
t940 = sin(t972);
t923 = -t940 * t978 + t941 * t982;
t922 = -t940 * t977 + t941 * t981;
t921 = -t940 * t976 + t941 * t980;
t920 = -t940 * t975 + t941 * t979;
t919 = -t940 * t982 - t941 * t978;
t918 = -t940 * t981 - t941 * t977;
t917 = -t940 * t980 - t941 * t976;
t916 = -t940 * t979 - t941 * t975;
t915 = t964 * t935 + t939 * t970;
t914 = -t935 * t970 + t939 * t964;
t913 = t962 * t934 + t938 * t968;
t912 = -t934 * t968 + t938 * t962;
t911 = t960 * t933 + t937 * t966;
t910 = -t933 * t966 + t937 * t960;
t909 = t952 * t932 + t936 * t954;
t908 = -t932 * t954 + t936 * t952;
t1 = [-m(4) * g(1) + (-t936 * t991 - t937 * t990 - t938 * t989 - t939 * t988) * t983 + (-t909 * t999 - t911 * t998 - t913 * t997 - t915 * t996) * t944; -m(4) * g(2) + (t932 * t991 + t933 * t990 + t934 * t989 + t935 * t988) * t983 + (-t908 * t999 - t910 * t998 - t912 * t997 - t914 * t996) * t944; -m(4) * g(3) + (-t943 * t924 * t995 - t948 * t925 * t994 - t949 * t926 * t993 - t950 * t927 * t992) * t944 + (-t970 / t969 ^ 2 * t907 * t992 - t968 / t967 ^ 2 * t906 * t993 - t966 / t965 ^ 2 * t905 * t994 - t954 / t953 ^ 2 * t904 * t995 + (t950 * ((t984 * rSges(3,2) - t1008) * t969 + t963 * (t984 * rSges(3,1) + t971)) + t949 * ((t985 * rSges(3,2) - t1008) * t967 + t961 * (t985 * rSges(3,1) + t971)) + t948 * ((t986 * rSges(3,2) - t1008) * t965 + t959 * (t986 * rSges(3,1) + t971)) + t943 * ((t987 * rSges(3,2) - t1008) * t953 + t951 * (t987 * rSges(3,1) + t971))) * m(3)) * t983; m(4) * ((g(1) * t974 + g(2) * t973) * t940 + (g(1) * t973 - g(2) * t974) * t941) + (-(t914 * t923 + t915 * t919) * t996 - (t912 * t922 + t913 * t918) * t997 - (t910 * t921 + t911 * t917) * t998 - (t908 * t920 + t909 * t916) * t999) * t944 + ((-t919 * t939 + t923 * t935) * t988 + (-t918 * t938 + t922 * t934) * t989 + (-t917 * t937 + t921 * t933) * t990 + (-t916 * t936 + t920 * t932) * t991) * t983;];
taugX  = t1;
