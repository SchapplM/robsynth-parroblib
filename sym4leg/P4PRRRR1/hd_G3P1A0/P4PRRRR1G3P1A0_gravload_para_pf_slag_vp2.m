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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
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

function taugX = P4PRRRR1G3P1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3P1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:03:07
% EndTime: 2020-03-02 19:03:08
% DurationCPUTime: 0.96s
% Computational Cost: add. (435->109), mult. (737->199), div. (92->13), fcn. (634->26), ass. (0->103)
t942 = mrSges(3,1) * g(3);
t893 = legFrame(4,2);
t869 = sin(t893);
t873 = cos(t893);
t861 = t869 * g(1) + t873 * g(2);
t889 = sin(qJ(2,4));
t879 = 0.1e1 / t889;
t941 = t861 * t879;
t894 = legFrame(3,2);
t870 = sin(t894);
t874 = cos(t894);
t862 = t870 * g(1) + t874 * g(2);
t898 = sin(qJ(2,3));
t882 = 0.1e1 / t898;
t940 = t862 * t882;
t895 = legFrame(2,2);
t871 = sin(t895);
t875 = cos(t895);
t863 = t871 * g(1) + t875 * g(2);
t900 = sin(qJ(2,2));
t883 = 0.1e1 / t900;
t939 = t863 * t883;
t896 = legFrame(1,2);
t872 = sin(t896);
t876 = cos(t896);
t864 = t872 * g(1) + t876 * g(2);
t902 = sin(qJ(2,1));
t884 = 0.1e1 / t902;
t938 = t864 * t884;
t888 = sin(qJ(3,4));
t937 = t879 * t888;
t897 = sin(qJ(3,3));
t936 = t882 * t897;
t899 = sin(qJ(3,2));
t935 = t883 * t899;
t901 = sin(qJ(3,1));
t934 = t884 * t901;
t865 = t873 * g(1) - t869 * g(2);
t891 = cos(qJ(2,4));
t892 = mrSges(2,2) - mrSges(3,3);
t890 = cos(qJ(3,4));
t925 = mrSges(3,1) * t890 - mrSges(3,2) * t888 + mrSges(2,1);
t841 = (t889 * t925 + t892 * t891) * t865 + (t889 * t892 - t925 * t891) * t861;
t880 = 0.1e1 / t890;
t933 = t841 * t879 * t880;
t866 = t874 * g(1) - t870 * g(2);
t904 = cos(qJ(2,3));
t903 = cos(qJ(3,3));
t924 = mrSges(3,1) * t903 - mrSges(3,2) * t897 + mrSges(2,1);
t842 = (t892 * t904 + t898 * t924) * t866 + (t898 * t892 - t924 * t904) * t862;
t885 = 0.1e1 / t903;
t932 = t842 * t882 * t885;
t867 = t875 * g(1) - t871 * g(2);
t906 = cos(qJ(2,2));
t905 = cos(qJ(3,2));
t923 = mrSges(3,1) * t905 - mrSges(3,2) * t899 + mrSges(2,1);
t843 = (t892 * t906 + t900 * t923) * t867 + (t900 * t892 - t923 * t906) * t863;
t886 = 0.1e1 / t905;
t931 = t843 * t883 * t886;
t868 = t876 * g(1) - t872 * g(2);
t908 = cos(qJ(2,1));
t907 = cos(qJ(3,1));
t922 = mrSges(3,1) * t907 - mrSges(3,2) * t901 + mrSges(2,1);
t844 = (t892 * t908 + t902 * t922) * t868 + (t902 * t892 - t922 * t908) * t864;
t887 = 0.1e1 / t907;
t930 = t844 * t884 * t887;
t929 = t861 * t889 + t865 * t891;
t928 = t862 * t898 + t866 * t904;
t927 = t863 * t900 + t867 * t906;
t926 = t864 * t902 + t868 * t908;
t921 = 0.1e1 / pkin(2);
t920 = koppelP(1,1);
t919 = koppelP(2,1);
t918 = koppelP(3,1);
t917 = koppelP(4,1);
t916 = koppelP(1,2);
t915 = koppelP(2,2);
t914 = koppelP(3,2);
t913 = koppelP(4,2);
t912 = mrSges(4,1);
t911 = mrSges(4,2);
t910 = xP(4);
t909 = mrSges(3,2) * g(3);
t881 = m(1) + m(2) + m(3);
t878 = cos(t910);
t877 = sin(t910);
t860 = -t877 * t916 + t878 * t920;
t859 = -t877 * t915 + t878 * t919;
t858 = -t877 * t914 + t878 * t918;
t857 = -t877 * t913 + t878 * t917;
t856 = -t877 * t920 - t878 * t916;
t855 = -t877 * t919 - t878 * t915;
t854 = -t877 * t918 - t878 * t914;
t853 = -t877 * t917 - t878 * t913;
t852 = t902 * t872 + t876 * t908;
t851 = -t872 * t908 + t876 * t902;
t850 = t900 * t871 + t875 * t906;
t849 = -t871 * t906 + t875 * t900;
t848 = t898 * t870 + t874 * t904;
t847 = -t870 * t904 + t874 * t898;
t846 = t889 * t869 + t873 * t891;
t845 = -t869 * t891 + t873 * t889;
t1 = [-g(1) * m(4) + (-t873 * t933 - t874 * t932 - t875 * t931 - t876 * t930) * t921 + (-t846 * t941 - t848 * t940 - t850 * t939 - t852 * t938) * t881; -g(2) * m(4) + (t869 * t933 + t870 * t932 + t871 * t931 + t872 * t930) * t921 + (-t845 * t941 - t847 * t940 - t849 * t939 - t851 * t938) * t881; -g(3) * m(4) + (-t880 * t861 * t937 - t885 * t862 * t936 - t886 * t863 * t935 - t887 * t864 * t934) * t881 + (-t908 / t907 ^ 2 * t844 * t934 + t887 * ((t926 * mrSges(3,2) - t942) * t907 + t901 * (t926 * mrSges(3,1) + t909)) - t906 / t905 ^ 2 * t843 * t935 + t886 * ((t927 * mrSges(3,2) - t942) * t905 + t899 * (t927 * mrSges(3,1) + t909)) - t904 / t903 ^ 2 * t842 * t936 + t885 * ((t928 * mrSges(3,2) - t942) * t903 + t897 * (t928 * mrSges(3,1) + t909)) - t891 / t890 ^ 2 * t841 * t937 + t880 * ((t929 * mrSges(3,2) - t942) * t890 + t888 * (t929 * mrSges(3,1) + t909))) * t921; -(-g(1) * t912 - g(2) * t911) * t877 + t878 * (g(1) * t911 - g(2) * t912) + (-(t851 * t860 + t852 * t856) * t938 - (t849 * t859 + t850 * t855) * t939 - (t847 * t858 + t848 * t854) * t940 - (t845 * t857 + t846 * t853) * t941) * t881 + ((-t856 * t876 + t860 * t872) * t930 + (-t855 * t875 + t859 * t871) * t931 + (-t854 * t874 + t858 * t870) * t932 + (-t853 * t873 + t857 * t869) * t933) * t921;];
taugX  = t1;
