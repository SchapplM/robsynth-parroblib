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
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V1G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:01:51
% EndTime: 2020-08-06 19:01:51
% DurationCPUTime: 0.35s
% Computational Cost: add. (525->110), mult. (735->185), div. (30->6), fcn. (540->18), ass. (0->93)
t864 = m(3) * pkin(1) + mrSges(2,1) + mrSges(3,1);
t919 = g(3) * t864;
t879 = sin(qJ(1,3));
t918 = t879 * pkin(4);
t881 = sin(qJ(1,2));
t917 = t881 * pkin(4);
t883 = sin(qJ(1,1));
t916 = t883 * pkin(4);
t915 = mrSges(3,3) - mrSges(2,2);
t878 = sin(qJ(2,3));
t884 = cos(qJ(2,3));
t891 = 0.1e1 / qJ(3,3);
t875 = legFrame(3,3);
t865 = sin(t875);
t868 = cos(t875);
t843 = -t865 * g(1) + t868 * g(2);
t846 = t868 * g(1) + t865 * g(2);
t885 = cos(qJ(1,3));
t899 = t843 * t879 + t846 * t885;
t914 = (-t884 * g(3) + t899 * t878) * t891;
t880 = sin(qJ(2,2));
t886 = cos(qJ(2,2));
t892 = 0.1e1 / qJ(3,2);
t876 = legFrame(2,3);
t866 = sin(t876);
t869 = cos(t876);
t844 = -t866 * g(1) + t869 * g(2);
t847 = t869 * g(1) + t866 * g(2);
t887 = cos(qJ(1,2));
t898 = t844 * t881 + t847 * t887;
t913 = (-t886 * g(3) + t898 * t880) * t892;
t882 = sin(qJ(2,1));
t888 = cos(qJ(2,1));
t893 = 0.1e1 / qJ(3,1);
t877 = legFrame(1,3);
t867 = sin(t877);
t870 = cos(t877);
t845 = -t867 * g(1) + t870 * g(2);
t848 = t870 * g(1) + t867 * g(2);
t889 = cos(qJ(1,1));
t897 = t845 * t883 + t848 * t889;
t912 = (-t888 * g(3) + t897 * t882) * t893;
t911 = t878 * qJ(3,3);
t910 = t880 * qJ(3,2);
t909 = t882 * qJ(3,1);
t890 = pkin(1) + pkin(2);
t908 = t890 * t884;
t907 = t890 * t886;
t906 = t890 * t888;
t861 = m(3) * qJ(3,3) + t915;
t905 = t891 * ((-t899 * t861 - t919) * t884 + (-g(3) * t861 + t899 * t864) * t878);
t862 = m(3) * qJ(3,2) + t915;
t904 = t892 * ((-t898 * t862 - t919) * t886 + (-g(3) * t862 + t898 * t864) * t880);
t863 = m(3) * qJ(3,1) + t915;
t903 = t893 * ((-t897 * t863 - t919) * t888 + (-g(3) * t863 + t897 * t864) * t882);
t902 = t884 * t905;
t901 = t886 * t904;
t900 = t888 * t903;
t896 = t861 * t878 + t864 * t884 + mrSges(1,1);
t895 = t862 * t880 + t864 * t886 + mrSges(1,1);
t894 = t863 * t882 + t864 * t888 + mrSges(1,1);
t874 = mrSges(1,2) - mrSges(3,2) - mrSges(2,3);
t873 = t889 * pkin(4);
t872 = t887 * pkin(4);
t871 = t885 * pkin(4);
t860 = t906 + t909;
t859 = t907 + t910;
t858 = t908 + t911;
t857 = 0.1e1 / t860;
t856 = 0.1e1 / t859;
t855 = 0.1e1 / t858;
t854 = t889 * t909 - t916;
t853 = t883 * t909 + t873;
t852 = t887 * t910 - t917;
t851 = t881 * t910 + t872;
t850 = t885 * t911 - t918;
t849 = t879 * t911 + t871;
t842 = t867 * t889 + t870 * t883;
t841 = -t867 * t883 + t870 * t889;
t840 = t866 * t887 + t869 * t881;
t839 = -t866 * t881 + t869 * t887;
t838 = t865 * t885 + t868 * t879;
t837 = -t865 * t879 + t868 * t885;
t836 = t860 * t889 - t916;
t835 = t859 * t887 - t917;
t834 = t858 * t885 - t918;
t833 = t860 * t883 + t873;
t832 = t859 * t881 + t872;
t831 = t858 * t879 + t871;
t824 = (t874 * t889 + t894 * t883) * t848 + (t874 * t883 - t894 * t889) * t845;
t823 = (t874 * t887 + t895 * t881) * t847 + (t874 * t881 - t895 * t887) * t844;
t822 = (t874 * t885 + t896 * t879) * t846 + (t874 * t879 - t896 * t885) * t843;
t1 = [-g(1) * m(4) + (-t842 * t824 + (t841 * t906 - t867 * t853 + t854 * t870) * t900) * t857 + (-t840 * t823 + (t839 * t907 - t866 * t851 + t852 * t869) * t901) * t856 + (-t838 * t822 + (t837 * t908 - t865 * t849 + t850 * t868) * t902) * t855 + (-(-t867 * t833 + t870 * t836) * t912 - (-t866 * t832 + t869 * t835) * t913 - (-t865 * t831 + t868 * t834) * t914) * m(3); -g(2) * m(4) + (t841 * t824 + (t842 * t906 + t853 * t870 + t867 * t854) * t900) * t857 + (t839 * t823 + (t840 * t907 + t851 * t869 + t866 * t852) * t901) * t856 + (t837 * t822 + (t838 * t908 + t849 * t868 + t865 * t850) * t902) * t855 + (-(t833 * t870 + t867 * t836) * t912 - (t832 * t869 + t866 * t835) * t913 - (t831 * t868 + t865 * t834) * t914) * m(3); t878 * t905 + t880 * t904 + t882 * t903 - g(3) * m(4) + (-(-t888 * qJ(3,1) + t890 * t882) * t912 - (-t886 * qJ(3,2) + t890 * t880) * t913 - (-t884 * qJ(3,3) + t890 * t878) * t914) * m(3);];
taugX  = t1;
