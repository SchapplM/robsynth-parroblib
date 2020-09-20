% Calculate Gravitation load for parallel robot
% P3RPRRR12V1G2A0
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR12V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:30
% EndTime: 2020-08-06 18:24:30
% DurationCPUTime: 0.44s
% Computational Cost: add. (423->112), mult. (645->193), div. (36->7), fcn. (375->18), ass. (0->77)
t754 = pkin(1) + pkin(5);
t753 = m(3) / pkin(3);
t717 = cos(qJ(3,3));
t752 = pkin(3) * t717;
t718 = cos(qJ(1,3));
t751 = pkin(3) * t718;
t719 = cos(qJ(3,2));
t750 = pkin(3) * t719;
t720 = cos(qJ(1,2));
t749 = pkin(3) * t720;
t721 = cos(qJ(3,1));
t748 = pkin(3) * t721;
t722 = cos(qJ(1,1));
t747 = pkin(3) * t722;
t708 = legFrame(3,2);
t698 = sin(t708);
t701 = cos(t708);
t683 = t701 * g(1) - t698 * g(2);
t712 = sin(qJ(1,3));
t675 = g(3) * t712 - t683 * t718;
t680 = t698 * g(1) + t701 * g(2);
t711 = sin(qJ(3,3));
t705 = 0.1e1 / t711;
t746 = ((t675 * rSges(3,1) - rSges(3,2) * t680) * t717 - t711 * (rSges(3,1) * t680 + t675 * rSges(3,2))) * t705;
t709 = legFrame(2,2);
t699 = sin(t709);
t702 = cos(t709);
t684 = t702 * g(1) - t699 * g(2);
t714 = sin(qJ(1,2));
t676 = g(3) * t714 - t684 * t720;
t681 = t699 * g(1) + t702 * g(2);
t713 = sin(qJ(3,2));
t706 = 0.1e1 / t713;
t745 = ((t676 * rSges(3,1) - rSges(3,2) * t681) * t719 - t713 * (rSges(3,1) * t681 + t676 * rSges(3,2))) * t706;
t710 = legFrame(1,2);
t700 = sin(t710);
t703 = cos(t710);
t685 = t703 * g(1) - t700 * g(2);
t716 = sin(qJ(1,1));
t677 = g(3) * t716 - t685 * t722;
t682 = t700 * g(1) + t703 * g(2);
t715 = sin(qJ(3,1));
t707 = 0.1e1 / t715;
t744 = ((t677 * rSges(3,1) - rSges(3,2) * t682) * t721 - t715 * (rSges(3,1) * t682 + t677 * rSges(3,2))) * t707;
t695 = t711 * pkin(3) + qJ(2,3);
t692 = 0.1e1 / t695;
t743 = t675 * t692;
t696 = t713 * pkin(3) + qJ(2,2);
t693 = 0.1e1 / t696;
t742 = t676 * t693;
t697 = t715 * pkin(3) + qJ(2,1);
t694 = 0.1e1 / t697;
t741 = t677 * t694;
t679 = (rSges(3,3) + t754) * m(3) + (pkin(1) - rSges(2,2)) * m(2) + m(1) * rSges(1,1);
t678 = t679 * g(3);
t723 = m(1) * rSges(1,2);
t689 = -qJ(2,3) * m(3) + (-rSges(2,3) - qJ(2,3)) * m(2) + t723;
t740 = t692 * ((t689 * g(3) - t683 * t679) * t718 + (t683 * t689 + t678) * t712 + (g(3) * t718 + t683 * t712) * m(3) * (-t711 * rSges(3,1) - rSges(3,2) * t717));
t690 = -qJ(2,2) * m(3) + (-rSges(2,3) - qJ(2,2)) * m(2) + t723;
t739 = t693 * ((t690 * g(3) - t684 * t679) * t720 + (t684 * t690 + t678) * t714 + (g(3) * t720 + t684 * t714) * m(3) * (-t713 * rSges(3,1) - rSges(3,2) * t719));
t691 = -qJ(2,1) * m(3) + (-rSges(2,3) - qJ(2,1)) * m(2) + t723;
t738 = t694 * ((t691 * g(3) - t685 * t679) * t722 + (t685 * t691 + t678) * t716 + (g(3) * t722 + t685 * t716) * m(3) * (-t715 * rSges(3,1) - rSges(3,2) * t721));
t737 = t717 * qJ(2,3);
t736 = t719 * qJ(2,2);
t735 = t721 * qJ(2,1);
t734 = t705 * t743;
t733 = t706 * t742;
t732 = t707 * t741;
t731 = t712 * t740;
t730 = t714 * t739;
t729 = t716 * t738;
t724 = -m(2) - m(3);
t704 = pkin(6) + t754;
t688 = qJ(2,1) * t722 - t704 * t716;
t687 = qJ(2,2) * t720 - t704 * t714;
t686 = qJ(2,3) * t718 - t704 * t712;
t1 = [t701 * t731 + t702 * t730 + t703 * t729 - m(4) * g(1) + (((-t688 * t703 + t700 * t748) * t715 + (t721 - 0.1e1) * (t721 + 0.1e1) * t703 * t747 + t700 * t735) * t732 + ((-t687 * t702 + t699 * t750) * t713 + (t719 - 0.1e1) * (t719 + 0.1e1) * t702 * t749 + t699 * t736) * t733 + ((-t686 * t701 + t698 * t752) * t711 + (t717 - 0.1e1) * (t717 + 0.1e1) * t701 * t751 + t698 * t737) * t734) * t724 + (t698 * t746 + t699 * t745 + t700 * t744) * t753; -t698 * t731 - t699 * t730 - t700 * t729 - m(4) * g(2) + (((t688 * t700 + t703 * t748) * t715 + (-t721 ^ 2 + 0.1e1) * t700 * t747 + t703 * t735) * t732 + ((t687 * t699 + t702 * t750) * t713 + (-t719 ^ 2 + 0.1e1) * t699 * t749 + t702 * t736) * t733 + ((t686 * t698 + t701 * t752) * t711 + (-t717 ^ 2 + 0.1e1) * t698 * t751 + t701 * t737) * t734) * t724 + (t701 * t746 + t702 * t745 + t703 * t744) * t753; t718 * t740 + t720 * t739 + t722 * t738 - m(4) * g(3) + ((t716 * t697 + t704 * t722) * t741 + (t714 * t696 + t704 * t720) * t742 + (t712 * t695 + t704 * t718) * t743) * t724;];
taugX  = t1;
