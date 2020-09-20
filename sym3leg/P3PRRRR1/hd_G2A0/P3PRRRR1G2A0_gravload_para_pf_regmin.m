% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR1G2A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:37
% EndTime: 2020-03-09 21:16:38
% DurationCPUTime: 0.46s
% Computational Cost: add. (180->63), mult. (426->153), div. (144->10), fcn. (504->18), ass. (0->77)
t926 = legFrame(3,2);
t905 = sin(t926);
t908 = cos(t926);
t899 = t905 * g(1) + t908 * g(2);
t930 = sin(qJ(2,3));
t936 = cos(qJ(2,3));
t890 = -g(3) * t930 + t899 * t936;
t915 = 0.1e1 / t930;
t935 = cos(qJ(3,3));
t920 = 0.1e1 / t935;
t974 = t915 * t920;
t986 = t890 * t974;
t927 = legFrame(2,2);
t906 = sin(t927);
t909 = cos(t927);
t900 = t906 * g(1) + t909 * g(2);
t932 = sin(qJ(2,2));
t938 = cos(qJ(2,2));
t891 = -g(3) * t932 + t900 * t938;
t917 = 0.1e1 / t932;
t937 = cos(qJ(3,2));
t922 = 0.1e1 / t937;
t972 = t917 * t922;
t985 = t891 * t972;
t928 = legFrame(1,2);
t907 = sin(t928);
t910 = cos(t928);
t901 = t907 * g(1) + t910 * g(2);
t934 = sin(qJ(2,1));
t940 = cos(qJ(2,1));
t892 = -g(3) * t934 + t901 * t940;
t919 = 0.1e1 / t934;
t939 = cos(qJ(3,1));
t924 = 0.1e1 / t939;
t970 = t919 * t924;
t984 = t892 * t970;
t902 = t908 * g(1) - t905 * g(2);
t929 = sin(qJ(3,3));
t959 = g(3) * t936 + t899 * t930;
t973 = t915 * t936;
t962 = t929 * t973;
t983 = t920 * (t890 * t962 + t902 * t935 + t929 * t959);
t903 = t909 * g(1) - t906 * g(2);
t931 = sin(qJ(3,2));
t958 = g(3) * t938 + t900 * t932;
t971 = t917 * t938;
t961 = t931 * t971;
t982 = t922 * (t891 * t961 + t903 * t937 + t931 * t958);
t904 = t910 * g(1) - t907 * g(2);
t933 = sin(qJ(3,1));
t957 = g(3) * t940 + t901 * t934;
t969 = t919 * t940;
t960 = t933 * t969;
t981 = t924 * (t892 * t960 + t904 * t939 + t933 * t957);
t968 = t930 * t935;
t967 = t932 * t937;
t966 = t934 * t939;
t965 = t899 * t974;
t964 = t900 * t972;
t963 = t901 * t970;
t921 = 0.1e1 / t935 ^ 2;
t956 = t921 * t962;
t923 = 0.1e1 / t937 ^ 2;
t955 = t923 * t961;
t925 = 0.1e1 / t939 ^ 2;
t954 = t925 * t960;
t953 = t905 * t956;
t952 = t906 * t955;
t951 = t907 * t954;
t950 = t908 * t956;
t949 = t909 * t955;
t948 = t910 * t954;
t944 = t890 * t929 ^ 2 * t921 * t973 - (-t902 * t929 + t935 * t959) * t920;
t943 = t891 * t931 ^ 2 * t923 * t971 - (-t903 * t931 + t937 * t958) * t922;
t942 = t892 * t933 ^ 2 * t925 * t969 - (-t904 * t933 + t939 * t957) * t924;
t941 = 0.1e1 / pkin(2);
t1 = [-(t907 * t966 - t910 * t933) * t963 - (t906 * t967 - t909 * t931) * t964 - (t905 * t968 - t908 * t929) * t965, 0, (-t890 * t950 - t891 * t949 - t892 * t948) * t941, (t948 * t957 + t949 * t958 + t950 * t959) * t941, 0, 0, 0, 0, 0, (-t908 * t983 - t909 * t982 - t910 * t981) * t941, (t944 * t908 + t943 * t909 + t942 * t910) * t941, -g(1); -(t907 * t933 + t910 * t966) * t963 - (t906 * t931 + t909 * t967) * t964 - (t905 * t929 + t908 * t968) * t965, 0, (t890 * t953 + t891 * t952 + t892 * t951) * t941, (-t951 * t957 - t952 * t958 - t953 * t959) * t941, 0, 0, 0, 0, 0, (t905 * t983 + t906 * t982 + t907 * t981) * t941, (-t944 * t905 - t943 * t906 - t942 * t907) * t941, -g(2); -t899 * t973 - t900 * t971 - t901 * t969, 0, (t984 + t985 + t986) * t941, (-t957 * t970 - t958 * t972 - t959 * t974) * t941, 0, 0, 0, 0, 0, (t890 * t915 + t891 * t917 + t892 * t919) * t941, (-t929 * t986 - t931 * t985 - t933 * t984) * t941, -g(3);];
tau_reg  = t1;
