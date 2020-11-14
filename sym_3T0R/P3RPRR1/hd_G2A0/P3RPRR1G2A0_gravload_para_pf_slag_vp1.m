% Calculate Gravitation load for parallel robot
% P3RPRR1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRR1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:24:51
% EndTime: 2020-03-09 21:24:51
% DurationCPUTime: 0.40s
% Computational Cost: add. (579->132), mult. (663->140), div. (27->4), fcn. (330->66), ass. (0->89)
t721 = m(2) + m(3);
t750 = m(1) * rSges(1,2);
t749 = m(2) * rSges(2,2);
t748 = m(3) / pkin(3);
t712 = legFrame(3,2);
t701 = sin(t712);
t704 = cos(t712);
t649 = t704 * g(1) - t701 * g(2);
t730 = pkin(7) + qJ(3,3);
t692 = qJ(1,3) + t730;
t685 = cos(t692);
t718 = rSges(3,2) * g(3);
t719 = rSges(3,1) * g(3);
t637 = -t685 * (rSges(3,1) * t649 - t718) + sin(t692) * (rSges(3,2) * t649 + t719);
t652 = 0.1e1 / (pkin(1) * sin(t730) + sin(qJ(3,3)) * pkin(2));
t747 = t637 * t652;
t713 = legFrame(2,2);
t702 = sin(t713);
t705 = cos(t713);
t650 = t705 * g(1) - t702 * g(2);
t731 = pkin(7) + qJ(3,2);
t693 = qJ(1,2) + t731;
t686 = cos(t693);
t638 = -t686 * (rSges(3,1) * t650 - t718) + sin(t693) * (rSges(3,2) * t650 + t719);
t653 = 0.1e1 / (pkin(1) * sin(t731) + sin(qJ(3,2)) * pkin(2));
t746 = t638 * t653;
t714 = legFrame(1,2);
t703 = sin(t714);
t706 = cos(t714);
t651 = t706 * g(1) - t703 * g(2);
t732 = pkin(7) + qJ(3,1);
t694 = qJ(1,1) + t732;
t687 = cos(t694);
t639 = -t687 * (rSges(3,1) * t651 - t718) + sin(t694) * (rSges(3,2) * t651 + t719);
t654 = 0.1e1 / (pkin(1) * sin(t732) + sin(qJ(3,1)) * pkin(2));
t745 = t639 * t654;
t744 = (t701 * g(1) + t704 * g(2)) * t721;
t743 = (t702 * g(1) + t705 * g(2)) * t721;
t742 = (t703 * g(1) + t706 * g(2)) * t721;
t668 = m(1) * rSges(1,1) + t721 * pkin(1);
t655 = t668 * g(3);
t688 = m(2) * rSges(2,1) + m(3) * pkin(2);
t675 = t688 * g(3);
t707 = qJ(1,3) + pkin(7);
t689 = cos(t707);
t710 = g(3) * t749;
t711 = g(3) * t750;
t715 = cos(qJ(1,3));
t741 = t652 * ((-t649 * t688 + t710) * t689 + (t649 * t749 + t675) * sin(t707) + (-t649 * t668 + t711) * t715 + (t649 * t750 + t655) * sin(qJ(1,3)) + t637 * m(3));
t708 = qJ(1,2) + pkin(7);
t690 = cos(t708);
t716 = cos(qJ(1,2));
t740 = t653 * ((-t650 * t688 + t710) * t690 + (t650 * t749 + t675) * sin(t708) + (-t650 * t668 + t711) * t716 + (t650 * t750 + t655) * sin(qJ(1,2)) + t638 * m(3));
t709 = qJ(1,1) + pkin(7);
t691 = cos(t709);
t717 = cos(qJ(1,1));
t739 = t654 * ((-t651 * t688 + t710) * t691 + (t651 * t749 + t675) * sin(t709) + (-t651 * t668 + t711) * t717 + (t651 * t750 + t655) * sin(qJ(1,1)) + t639 * m(3));
t676 = t712 + t707;
t669 = qJ(3,3) + t676;
t677 = -t712 + t707;
t670 = qJ(3,3) + t677;
t738 = sin(t669) + sin(t670);
t678 = t713 + t708;
t671 = qJ(3,2) + t678;
t679 = -t713 + t708;
t672 = qJ(3,2) + t679;
t737 = sin(t671) + sin(t672);
t680 = t714 + t709;
t673 = qJ(3,1) + t680;
t681 = -t714 + t709;
t674 = qJ(3,1) + t681;
t736 = sin(t673) + sin(t674);
t735 = -cos(t670) + cos(t669);
t734 = -cos(t672) + cos(t671);
t733 = -cos(t674) + cos(t673);
t729 = t748 / 0.2e1;
t728 = t741 / 0.2e1;
t727 = t740 / 0.2e1;
t726 = t739 / 0.2e1;
t725 = t729 * t747;
t724 = t729 * t746;
t723 = t729 * t745;
t700 = qJ(1,1) - t714;
t699 = qJ(1,1) + t714;
t698 = qJ(1,2) - t713;
t697 = qJ(1,2) + t713;
t696 = qJ(1,3) - t712;
t695 = qJ(1,3) + t712;
t1 = [t736 * t726 - t703 * t742 - (t736 * pkin(3) + (sin(t680) + sin(t681)) * pkin(2) + (sin(t699) + sin(t700)) * pkin(1)) * t723 + t737 * t727 - t702 * t743 - (t737 * pkin(3) + (sin(t678) + sin(t679)) * pkin(2) + (sin(t697) + sin(t698)) * pkin(1)) * t724 + t738 * t728 - t701 * t744 - (t738 * pkin(3) + (sin(t676) + sin(t677)) * pkin(2) + (sin(t695) + sin(t696)) * pkin(1)) * t725 - m(4) * g(1); t733 * t726 - t706 * t742 + (-t733 * pkin(3) + (cos(t681) - cos(t680)) * pkin(2) + (cos(t700) - cos(t699)) * pkin(1)) * t723 + t734 * t727 - t705 * t743 + (-t734 * pkin(3) + (cos(t679) - cos(t678)) * pkin(2) + (cos(t698) - cos(t697)) * pkin(1)) * t724 + t735 * t728 - t704 * t744 + (-t735 * pkin(3) + (cos(t677) - cos(t676)) * pkin(2) + (cos(t696) - cos(t695)) * pkin(1)) * t725 - m(4) * g(2); t685 * t741 + t686 * t740 + t687 * t739 - m(4) * g(3) + ((-t717 * pkin(1) - pkin(2) * t691 - pkin(3) * t687) * t745 + (-t716 * pkin(1) - pkin(2) * t690 - pkin(3) * t686) * t746 + (-t715 * pkin(1) - pkin(2) * t689 - pkin(3) * t685) * t747) * t748;];
taugX  = t1;
