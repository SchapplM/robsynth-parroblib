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
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V1G3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:10:08
% EndTime: 2020-08-06 19:10:09
% DurationCPUTime: 0.79s
% Computational Cost: add. (591->119), mult. (786->209), div. (36->6), fcn. (534->18), ass. (0->96)
t858 = cos(qJ(2,3));
t865 = pkin(1) + pkin(2);
t852 = sin(qJ(2,3));
t886 = t852 * qJ(3,3);
t829 = t865 * t858 + t886;
t853 = sin(qJ(1,3));
t859 = cos(qJ(1,3));
t911 = t859 * pkin(4) + t829 * t853;
t905 = t853 * pkin(4);
t910 = t829 * t859 - t905;
t860 = cos(qJ(2,2));
t854 = sin(qJ(2,2));
t884 = t854 * qJ(3,2);
t830 = t865 * t860 + t884;
t855 = sin(qJ(1,2));
t861 = cos(qJ(1,2));
t909 = t861 * pkin(4) + t830 * t855;
t904 = t855 * pkin(4);
t908 = t830 * t861 - t904;
t862 = cos(qJ(2,1));
t856 = sin(qJ(2,1));
t882 = t856 * qJ(3,1);
t831 = t865 * t862 + t882;
t857 = sin(qJ(1,1));
t863 = cos(qJ(1,1));
t907 = t863 * pkin(4) + t831 * t857;
t903 = t857 * pkin(4);
t906 = t831 * t863 - t903;
t899 = mrSges(3,3) - mrSges(2,2);
t849 = legFrame(3,2);
t839 = sin(t849);
t842 = cos(t849);
t814 = t839 * g(1) + t842 * g(2);
t835 = m(3) * qJ(3,3) + t899;
t838 = m(3) * pkin(1) + mrSges(2,1) + mrSges(3,1);
t866 = 0.1e1 / qJ(3,3);
t817 = t842 * g(1) - t839 * g(2);
t874 = g(3) * t853 - t817 * t859;
t898 = ((-t814 * t838 + t874 * t835) * t858 + t852 * (-t814 * t835 - t874 * t838)) * t866;
t850 = legFrame(2,2);
t840 = sin(t850);
t843 = cos(t850);
t815 = t840 * g(1) + t843 * g(2);
t836 = m(3) * qJ(3,2) + t899;
t867 = 0.1e1 / qJ(3,2);
t818 = t843 * g(1) - t840 * g(2);
t873 = g(3) * t855 - t818 * t861;
t897 = ((-t815 * t838 + t873 * t836) * t860 + t854 * (-t815 * t836 - t873 * t838)) * t867;
t851 = legFrame(1,2);
t841 = sin(t851);
t844 = cos(t851);
t816 = t841 * g(1) + t844 * g(2);
t837 = m(3) * qJ(3,1) + t899;
t868 = 0.1e1 / qJ(3,1);
t819 = t844 * g(1) - t841 * g(2);
t872 = g(3) * t857 - t819 * t863;
t896 = ((-t816 * t838 + t872 * t837) * t862 + t856 * (-t816 * t837 - t872 * t838)) * t868;
t895 = (-t858 * t814 - t852 * t874) * t866;
t894 = (-t860 * t815 - t854 * t873) * t867;
t893 = (-t862 * t816 - t856 * t872) * t868;
t892 = t839 * qJ(3,3);
t891 = t840 * qJ(3,2);
t890 = t841 * qJ(3,1);
t889 = t842 * qJ(3,3);
t888 = t843 * qJ(3,2);
t887 = t844 * qJ(3,1);
t845 = mrSges(1,2) - mrSges(3,2) - mrSges(2,3);
t864 = g(3) * mrSges(1,1);
t871 = t852 * t835 + t838 * t858;
t805 = t864 * t859 + (-t845 * t853 + t871 * t859) * g(3) + (t845 * t859 + (mrSges(1,1) + t871) * t853) * t817;
t885 = t853 * t805;
t870 = t854 * t836 + t838 * t860;
t806 = t864 * t861 + (-t845 * t855 + t870 * t861) * g(3) + (t845 * t861 + (mrSges(1,1) + t870) * t855) * t818;
t883 = t855 * t806;
t869 = t856 * t837 + t838 * t862;
t807 = t864 * t863 + (-t845 * t857 + t869 * t863) * g(3) + (t845 * t863 + (mrSges(1,1) + t869) * t857) * t819;
t881 = t857 * t807;
t880 = t859 * t865;
t879 = t861 * t865;
t878 = t863 * t865;
t877 = t865 * t852;
t876 = t865 * t854;
t875 = t865 * t856;
t848 = t862 ^ 2;
t847 = t860 ^ 2;
t846 = t858 ^ 2;
t828 = -t862 * qJ(3,1) + t875;
t827 = -t860 * qJ(3,2) + t876;
t826 = -t858 * qJ(3,3) + t877;
t825 = 0.1e1 / t831;
t824 = 0.1e1 / t830;
t823 = 0.1e1 / t829;
t822 = t863 * t882 - t903;
t821 = t861 * t884 - t904;
t820 = t859 * t886 - t905;
t1 = [-g(1) * m(4) + (-t844 * t881 + ((t844 * t878 - t890) * t848 + (t822 * t844 + t841 * t875) * t862 + t890) * t896) * t825 + (-t843 * t883 + ((t843 * t879 - t891) * t847 + (t821 * t843 + t840 * t876) * t860 + t891) * t897) * t824 + (-t842 * t885 + ((t842 * t880 - t892) * t846 + (t820 * t842 + t839 * t877) * t858 + t892) * t898) * t823 + (-(t828 * t841 + t906 * t844) * t893 - (t827 * t840 + t908 * t843) * t894 - (t826 * t839 + t910 * t842) * t895) * m(3); -g(2) * m(4) + (t841 * t881 + ((-t841 * t878 - t887) * t848 + (-t841 * t822 + t844 * t875) * t862 + t887) * t896) * t825 + (t840 * t883 + ((-t840 * t879 - t888) * t847 + (-t840 * t821 + t843 * t876) * t860 + t888) * t897) * t824 + (t839 * t885 + ((-t839 * t880 - t889) * t846 + (-t839 * t820 + t842 * t877) * t858 + t889) * t898) * t823 + (-(t828 * t844 - t906 * t841) * t893 - (t827 * t843 - t908 * t840) * t894 - (t826 * t842 - t910 * t839) * t895) * m(3); -g(3) * m(4) + (-t907 * t862 * t896 - t863 * t807) * t825 + (-t909 * t860 * t897 - t861 * t806) * t824 + (-t911 * t858 * t898 - t859 * t805) * t823 + (t907 * t893 + t909 * t894 + t911 * t895) * m(3);];
taugX  = t1;
