% Calculate Gravitation load for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
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
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR1V1G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:08:00
% EndTime: 2022-11-04 17:08:00
% DurationCPUTime: 0.26s
% Computational Cost: add. (426->95), mult. (669->169), div. (36->7), fcn. (342->18), ass. (0->77)
t698 = m(1) * rSges(1,1);
t650 = legFrame(3,2);
t633 = sin(t650);
t636 = cos(t650);
t622 = t633 * g(1) + t636 * g(2);
t631 = (pkin(1) + rSges(3,1)) * m(3) + m(2) * rSges(2,1);
t632 = m(2) * rSges(2,2) + m(3) * rSges(3,2);
t659 = cos(qJ(2,3));
t643 = 0.1e1 / t659;
t653 = sin(qJ(2,3));
t625 = t636 * g(1) - t633 * g(2);
t654 = sin(qJ(1,3));
t660 = cos(qJ(1,3));
t675 = g(3) * t660 + t625 * t654;
t697 = ((t622 * t632 + t675 * t631) * t653 + t659 * (-t622 * t631 + t675 * t632)) * t643;
t651 = legFrame(2,2);
t634 = sin(t651);
t637 = cos(t651);
t623 = t634 * g(1) + t637 * g(2);
t661 = cos(qJ(2,2));
t644 = 0.1e1 / t661;
t655 = sin(qJ(2,2));
t626 = t637 * g(1) - t634 * g(2);
t656 = sin(qJ(1,2));
t662 = cos(qJ(1,2));
t674 = g(3) * t662 + t626 * t656;
t696 = ((t623 * t632 + t674 * t631) * t655 + t661 * (-t623 * t631 + t674 * t632)) * t644;
t652 = legFrame(1,2);
t635 = sin(t652);
t638 = cos(t652);
t624 = t635 * g(1) + t638 * g(2);
t663 = cos(qJ(2,1));
t645 = 0.1e1 / t663;
t657 = sin(qJ(2,1));
t627 = t638 * g(1) - t635 * g(2);
t658 = sin(qJ(1,1));
t664 = cos(qJ(1,1));
t673 = g(3) * t664 + t627 * t658;
t695 = ((t624 * t632 + t673 * t631) * t657 + t663 * (-t624 * t631 + t673 * t632)) * t645;
t647 = pkin(3) + qJ(3,3);
t640 = 0.1e1 / t647;
t694 = (g(3) * t654 - t625 * t660) * t640;
t648 = pkin(3) + qJ(3,2);
t641 = 0.1e1 / t648;
t693 = (g(3) * t656 - t626 * t662) * t641;
t649 = pkin(3) + qJ(3,1);
t642 = 0.1e1 / t649;
t692 = (g(3) * t658 - t627 * t664) * t642;
t676 = m(1) * rSges(1,2) - rSges(2,3) * m(2);
t628 = (-rSges(3,3) - qJ(3,3)) * m(3) + t676;
t639 = g(3) * t698;
t672 = t659 * t631 - t632 * t653;
t691 = t640 * (t639 * t654 + (t628 * t660 + t672 * t654) * g(3) + ((-t672 - t698) * t660 + t628 * t654) * t625);
t629 = (-rSges(3,3) - qJ(3,2)) * m(3) + t676;
t671 = t661 * t631 - t632 * t655;
t690 = t641 * (t639 * t656 + (t629 * t662 + t671 * t656) * g(3) + ((-t671 - t698) * t662 + t629 * t656) * t626);
t630 = (-rSges(3,3) - qJ(3,1)) * m(3) + t676;
t670 = t663 * t631 - t632 * t657;
t689 = t642 * (t639 * t658 + (t630 * t664 + t670 * t658) * g(3) + ((-t670 - t698) * t664 + t630 * t658) * t627);
t666 = pkin(2) + pkin(1);
t688 = t653 * t666;
t687 = t654 * t659;
t686 = t655 * t666;
t685 = t656 * t661;
t684 = t657 * t666;
t683 = t658 * t663;
t682 = t659 * t666;
t681 = t661 * t666;
t680 = t663 * t666;
t679 = t643 * t691;
t678 = t644 * t690;
t677 = t645 * t689;
t669 = -t647 * t660 + t654 * t682;
t668 = -t648 * t662 + t656 * t681;
t667 = -t649 * t664 + t658 * t680;
t646 = 0.1e1 / t666;
t1 = [(t635 * t657 + t638 * t683) * t677 + (t634 * t655 + t637 * t685) * t678 + (t633 * t653 + t636 * t687) * t679 - m(4) * g(1) + (t633 * t697 + t634 * t696 + t635 * t695) * t646 + (-(t635 * t684 + t667 * t638) * t692 - (t634 * t686 + t668 * t637) * t693 - (t633 * t688 + t669 * t636) * t694) * m(3); (-t635 * t683 + t657 * t638) * t677 + (-t634 * t685 + t655 * t637) * t678 + (-t633 * t687 + t653 * t636) * t679 - m(4) * g(2) + (t636 * t697 + t637 * t696 + t638 * t695) * t646 + (-(-t667 * t635 + t638 * t684) * t692 - (-t668 * t634 + t637 * t686) * t693 - (-t669 * t633 + t636 * t688) * t694) * m(3); t660 * t691 + t662 * t690 + t664 * t689 - m(4) * g(3) + (-(t658 * t649 + t664 * t680) * t692 - (t656 * t648 + t662 * t681) * t693 - (t654 * t647 + t660 * t682) * t694) * m(3);];
taugX  = t1;
