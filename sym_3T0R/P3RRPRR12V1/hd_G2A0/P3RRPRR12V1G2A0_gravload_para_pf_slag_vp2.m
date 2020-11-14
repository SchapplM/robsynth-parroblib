% Calculate Gravitation load for parallel robot
% P3RRPRR12V1G2A0
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
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:05:49
% EndTime: 2020-08-06 19:05:49
% DurationCPUTime: 0.49s
% Computational Cost: add. (591->119), mult. (792->207), div. (36->6), fcn. (540->18), ass. (0->97)
t888 = mrSges(3,3) - mrSges(2,2);
t819 = m(3) * qJ(3,1) + t888;
t821 = m(3) * pkin(1) + mrSges(2,1) + mrSges(3,1);
t839 = sin(qJ(2,1));
t845 = cos(qJ(2,1));
t894 = t839 * t819 + t821 * t845;
t818 = m(3) * qJ(3,2) + t888;
t837 = sin(qJ(2,2));
t843 = cos(qJ(2,2));
t893 = t837 * t818 + t821 * t843;
t817 = m(3) * qJ(3,3) + t888;
t835 = sin(qJ(2,3));
t841 = cos(qJ(2,3));
t892 = t835 * t817 + t821 * t841;
t842 = cos(qJ(1,3));
t891 = t842 * pkin(4);
t844 = cos(qJ(1,2));
t890 = t844 * pkin(4);
t846 = cos(qJ(1,1));
t889 = t846 * pkin(4);
t832 = legFrame(3,2);
t822 = sin(t832);
t825 = cos(t832);
t799 = t822 * g(1) + t825 * g(2);
t849 = 0.1e1 / qJ(3,3);
t802 = t825 * g(1) - t822 * g(2);
t836 = sin(qJ(1,3));
t854 = g(3) * t842 + t802 * t836;
t887 = ((-t799 * t821 - t854 * t817) * t841 + t835 * (-t799 * t817 + t854 * t821)) * t849;
t833 = legFrame(2,2);
t823 = sin(t833);
t826 = cos(t833);
t800 = t823 * g(1) + t826 * g(2);
t850 = 0.1e1 / qJ(3,2);
t803 = t826 * g(1) - t823 * g(2);
t838 = sin(qJ(1,2));
t853 = g(3) * t844 + t803 * t838;
t886 = ((-t800 * t821 - t853 * t818) * t843 + t837 * (-t800 * t818 + t853 * t821)) * t850;
t834 = legFrame(1,2);
t824 = sin(t834);
t827 = cos(t834);
t801 = t824 * g(1) + t827 * g(2);
t851 = 0.1e1 / qJ(3,1);
t804 = t827 * g(1) - t824 * g(2);
t840 = sin(qJ(1,1));
t852 = g(3) * t846 + t804 * t840;
t885 = ((-t801 * t821 - t852 * t819) * t845 + t839 * (-t801 * t819 + t852 * t821)) * t851;
t884 = (-t841 * t799 + t835 * t854) * t849;
t883 = (-t843 * t800 + t837 * t853) * t850;
t882 = (-t845 * t801 + t839 * t852) * t851;
t878 = t822 * qJ(3,3);
t877 = t823 * qJ(3,2);
t876 = t824 * qJ(3,1);
t875 = t825 * qJ(3,3);
t874 = t826 * qJ(3,2);
t873 = t827 * qJ(3,1);
t872 = t835 * qJ(3,3);
t848 = pkin(1) + pkin(2);
t870 = t836 * t848;
t869 = t837 * qJ(3,2);
t867 = t838 * t848;
t866 = t839 * qJ(3,1);
t864 = t840 * t848;
t828 = mrSges(1,2) - mrSges(3,2) - mrSges(2,3);
t820 = t828 * g(3);
t847 = g(3) * mrSges(1,1);
t787 = t820 * t842 + t836 * (t892 * g(3) + t847) + ((-mrSges(1,1) - t892) * t842 + t836 * t828) * t802;
t863 = t842 * t787;
t788 = t820 * t844 + t838 * (t893 * g(3) + t847) + ((-mrSges(1,1) - t893) * t844 + t838 * t828) * t803;
t862 = t844 * t788;
t789 = t820 * t846 + t840 * (t894 * g(3) + t847) + ((-mrSges(1,1) - t894) * t846 + t840 * t828) * t804;
t861 = t846 * t789;
t860 = t848 * t835;
t859 = t848 * t837;
t858 = t848 * t839;
t814 = t848 * t841 + t872;
t857 = -t836 * pkin(4) + t814 * t842;
t815 = t848 * t843 + t869;
t856 = -t838 * pkin(4) + t815 * t844;
t816 = t848 * t845 + t866;
t855 = -t840 * pkin(4) + t816 * t846;
t831 = t845 ^ 2;
t830 = t843 ^ 2;
t829 = t841 ^ 2;
t813 = -t845 * qJ(3,1) + t858;
t812 = -t843 * qJ(3,2) + t859;
t811 = -t841 * qJ(3,3) + t860;
t810 = 0.1e1 / t816;
t809 = 0.1e1 / t815;
t808 = 0.1e1 / t814;
t807 = t840 * t866 + t889;
t806 = t838 * t869 + t890;
t805 = t836 * t872 + t891;
t798 = t816 * t840 + t889;
t797 = t815 * t838 + t890;
t796 = t814 * t836 + t891;
t1 = [-g(1) * m(4) + (t827 * t861 + ((t827 * t864 - t876) * t831 + (t807 * t827 + t824 * t858) * t845 + t876) * t885) * t810 + (t826 * t862 + ((t826 * t867 - t877) * t830 + (t806 * t826 + t823 * t859) * t843 + t877) * t886) * t809 + (t825 * t863 + ((t825 * t870 - t878) * t829 + (t805 * t825 + t822 * t860) * t841 + t878) * t887) * t808 + (-(t798 * t827 + t824 * t813) * t882 - (t797 * t826 + t823 * t812) * t883 - (t796 * t825 + t822 * t811) * t884) * m(3); -g(2) * m(4) + (-t824 * t861 + ((-t824 * t864 - t873) * t831 + (-t824 * t807 + t827 * t858) * t845 + t873) * t885) * t810 + (-t823 * t862 + ((-t823 * t867 - t874) * t830 + (-t823 * t806 + t826 * t859) * t843 + t874) * t886) * t809 + (-t822 * t863 + ((-t822 * t870 - t875) * t829 + (-t822 * t805 + t825 * t860) * t841 + t875) * t887) * t808 + (-(-t798 * t824 + t813 * t827) * t882 - (-t797 * t823 + t812 * t826) * t883 - (-t796 * t822 + t811 * t825) * t884) * m(3); -g(3) * m(4) + (t845 * t855 * t885 - t840 * t789) * t810 + (t843 * t856 * t886 - t838 * t788) * t809 + (t841 * t857 * t887 - t836 * t787) * t808 + (-t855 * t882 - t856 * t883 - t857 * t884) * m(3);];
taugX  = t1;
