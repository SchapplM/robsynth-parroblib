% Calculate Gravitation load for parallel robot
% P3PRRRR8V1G3A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR8V1G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:16:00
% EndTime: 2020-08-06 17:16:01
% DurationCPUTime: 1.09s
% Computational Cost: add. (675->128), mult. (1692->253), div. (54->7), fcn. (1311->22), ass. (0->109)
t930 = legFrame(3,2);
t916 = sin(t930);
t919 = cos(t930);
t898 = t916 * g(1) + t919 * g(2);
t927 = sin(pkin(3));
t929 = cos(pkin(3));
t901 = t919 * g(1) - t916 * g(2);
t928 = cos(pkin(6));
t915 = g(3) * t928;
t926 = sin(pkin(6));
t967 = -t901 * t926 - t915;
t1020 = t898 * t927 + t967 * t929;
t934 = sin(qJ(2,3));
t940 = cos(qJ(2,3));
t1006 = g(3) * t926;
t966 = -t901 * t928 + t1006;
t949 = t1020 * t934 - t966 * t940;
t931 = legFrame(2,2);
t917 = sin(t931);
t920 = cos(t931);
t899 = t917 * g(1) + t920 * g(2);
t902 = t920 * g(1) - t917 * g(2);
t965 = -t902 * t926 - t915;
t1021 = t899 * t927 + t965 * t929;
t936 = sin(qJ(2,2));
t942 = cos(qJ(2,2));
t964 = -t902 * t928 + t1006;
t948 = t1021 * t936 - t964 * t942;
t932 = legFrame(1,2);
t918 = sin(t932);
t921 = cos(t932);
t900 = t918 * g(1) + t921 * g(2);
t903 = t921 * g(1) - t918 * g(2);
t963 = -t903 * t926 - t915;
t1022 = t900 * t927 + t963 * t929;
t938 = sin(qJ(2,1));
t944 = cos(qJ(2,1));
t962 = -t903 * t928 + t1006;
t947 = t1022 * t938 - t962 * t944;
t933 = sin(qJ(3,3));
t1012 = pkin(2) * t933;
t939 = cos(qJ(3,3));
t1009 = pkin(2) * t939;
t904 = -t940 * pkin(5) + t934 * t1009;
t889 = t929 * t1012 + t904 * t927;
t1016 = 0.1e1 / t889;
t935 = sin(qJ(3,2));
t1011 = pkin(2) * t935;
t941 = cos(qJ(3,2));
t1008 = pkin(2) * t941;
t905 = -t942 * pkin(5) + t936 * t1008;
t890 = t929 * t1011 + t905 * t927;
t1015 = 0.1e1 / t890;
t937 = sin(qJ(3,1));
t1010 = pkin(2) * t937;
t943 = cos(qJ(3,1));
t1007 = pkin(2) * t943;
t906 = -t944 * pkin(5) + t938 * t1007;
t891 = t929 * t1010 + t906 * t927;
t1014 = 0.1e1 / t891;
t1013 = m(3) / pkin(2);
t1005 = t1016 / t939;
t1004 = t1015 / t941;
t1003 = t1014 / t943;
t1002 = t1016 * t898;
t1001 = t1015 * t899;
t1000 = t1014 * t900;
t998 = t898 * t929;
t996 = t899 * t929;
t994 = t900 * t929;
t990 = rSges(3,2) * t927 * t915;
t989 = t926 * t927;
t988 = t929 * t934;
t987 = t929 * t936;
t986 = t929 * t938;
t985 = t929 * t940;
t984 = t929 * t942;
t983 = t929 * t944;
t982 = t933 * t940;
t981 = t935 * t942;
t980 = t937 * t944;
t979 = (((t967 * t927 - t998) * rSges(3,1) + t949 * rSges(3,2)) * t939 + t933 * (t990 + (t901 * t989 + t998) * rSges(3,2) + t949 * rSges(3,1))) * t1005;
t978 = (((t965 * t927 - t996) * rSges(3,1) + t948 * rSges(3,2)) * t941 + t935 * (t990 + (t902 * t989 + t996) * rSges(3,2) + t948 * rSges(3,1))) * t1004;
t977 = (((t963 * t927 - t994) * rSges(3,1) + t947 * rSges(3,2)) * t943 + t937 * (t990 + (t903 * t989 + t994) * rSges(3,2) + t947 * rSges(3,1))) * t1003;
t913 = m(2) * rSges(2,2) - m(3) * rSges(3,3);
t945 = m(2) * rSges(2,1);
t976 = (t949 * t913 + (-t1020 * t940 - t966 * t934) * (t945 + (rSges(3,1) * t939 - rSges(3,2) * t933) * m(3))) * t1005;
t975 = (t948 * t913 + (-t1021 * t942 - t964 * t936) * (t945 + (rSges(3,1) * t941 - rSges(3,2) * t935) * m(3))) * t1004;
t974 = (t947 * t913 + (-t1022 * t944 - t962 * t938) * (t945 + (rSges(3,1) * t943 - rSges(3,2) * t937) * m(3))) * t1003;
t973 = ((-t926 * t934 + t928 * t985) * t1009 + pkin(5) * (t926 * t940 + t928 * t988)) * t979;
t972 = ((-t926 * t936 + t928 * t984) * t1008 + pkin(5) * (t926 * t942 + t928 * t987)) * t978;
t971 = ((-t926 * t938 + t928 * t983) * t1007 + pkin(5) * (t926 * t944 + t928 * t986)) * t977;
t955 = t927 * t939 + t933 * t988;
t970 = (t926 * t982 + t955 * t928) * t976;
t954 = t927 * t941 + t935 * t987;
t969 = (t926 * t981 + t954 * t928) * t975;
t953 = t927 * t943 + t937 * t986;
t968 = (t926 * t980 + t953 * t928) * t974;
t958 = t927 * t1012 - t904 * t929;
t957 = t927 * t1011 - t905 * t929;
t956 = t927 * t1010 - t906 * t929;
t922 = m(1) + m(2) + m(3);
t909 = pkin(5) * t938 + t944 * t1007;
t908 = pkin(5) * t936 + t942 * t1008;
t907 = pkin(5) * t934 + t940 * t1009;
t876 = t928 * t909 + t956 * t926;
t875 = t928 * t908 + t957 * t926;
t874 = t928 * t907 + t958 * t926;
t1 = [-t919 * t970 - t920 * t969 - t921 * t968 - m(4) * g(1) + (-(t876 * t921 + t918 * t891) * t1000 - (t875 * t920 + t917 * t890) * t1001 - (t874 * t919 + t916 * t889) * t1002) * t922 + (-t919 * t973 - t920 * t972 - t921 * t971) * t1013; t916 * t970 + t917 * t969 + t918 * t968 - m(4) * g(2) + (-(-t876 * t918 + t921 * t891) * t1000 - (-t875 * t917 + t920 * t890) * t1001 - (-t874 * t916 + t919 * t889) * t1002) * t922 + (t916 * t973 + t917 * t972 + t918 * t971) * t1013; (t953 * t926 - t928 * t980) * t974 + (t954 * t926 - t928 * t981) * t975 + (t955 * t926 - t928 * t982) * t976 - m(4) * g(3) + (-(-t926 * t909 + t956 * t928) * t1000 - (-t926 * t908 + t957 * t928) * t1001 - (-t926 * t907 + t958 * t928) * t1002) * t922 + (((t926 * t983 + t928 * t938) * t1007 + (t926 * t986 - t928 * t944) * pkin(5)) * t977 + ((t926 * t984 + t928 * t936) * t1008 + (t926 * t987 - t928 * t942) * pkin(5)) * t978 + ((t926 * t985 + t928 * t934) * t1009 + (t926 * t988 - t928 * t940) * pkin(5)) * t979) * t1013;];
taugX  = t1;
