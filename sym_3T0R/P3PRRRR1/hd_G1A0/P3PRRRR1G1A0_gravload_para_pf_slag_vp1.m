% Calculate Gravitation load for parallel robot
% P3PRRRR1G1A0
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
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR1G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:03
% EndTime: 2020-03-09 20:34:04
% DurationCPUTime: 0.31s
% Computational Cost: add. (183->74), mult. (401->149), div. (42->10), fcn. (276->18), ass. (0->61)
t666 = sin(qJ(3,1));
t672 = cos(qJ(3,1));
t705 = -rSges(3,1) * t672 + rSges(3,2) * t666;
t664 = sin(qJ(3,2));
t670 = cos(qJ(3,2));
t704 = -rSges(3,1) * t670 + rSges(3,2) * t664;
t662 = sin(qJ(3,3));
t668 = cos(qJ(3,3));
t703 = -rSges(3,1) * t668 + rSges(3,2) * t662;
t702 = rSges(2,1) * g(3);
t701 = g(3) * rSges(3,3);
t649 = (m(1) + m(2) + m(3));
t700 = g(3) * t649;
t659 = legFrame(3,3);
t643 = sin(t659);
t646 = cos(t659);
t637 = -t643 * g(1) + t646 * g(2);
t653 = 0.1e1 / t668;
t640 = t646 * g(1) + t643 * g(2);
t663 = sin(qJ(2,3));
t669 = cos(qJ(2,3));
t678 = g(3) * t663 + t640 * t669;
t693 = ((rSges(3,1) * t637 + t678 * rSges(3,2)) * t668 + t662 * (t678 * rSges(3,1) - rSges(3,2) * t637)) * t653;
t660 = legFrame(2,3);
t644 = sin(t660);
t647 = cos(t660);
t638 = -t644 * g(1) + t647 * g(2);
t655 = 0.1e1 / t670;
t641 = t647 * g(1) + t644 * g(2);
t665 = sin(qJ(2,2));
t671 = cos(qJ(2,2));
t677 = g(3) * t665 + t641 * t671;
t692 = ((rSges(3,1) * t638 + t677 * rSges(3,2)) * t670 + t664 * (t677 * rSges(3,1) - rSges(3,2) * t638)) * t655;
t661 = legFrame(1,3);
t645 = sin(t661);
t648 = cos(t661);
t639 = -t645 * g(1) + t648 * g(2);
t657 = 0.1e1 / t672;
t642 = t648 * g(1) + t645 * g(2);
t667 = sin(qJ(2,1));
t673 = cos(qJ(2,1));
t676 = g(3) * t667 + t642 * t673;
t691 = ((rSges(3,1) * t639 + t676 * rSges(3,2)) * t672 + t666 * (t676 * rSges(3,1) - rSges(3,2) * t639)) * t657;
t650 = 0.1e1 / t663;
t690 = t650 * t653;
t651 = 0.1e1 / t665;
t689 = t651 * t655;
t652 = 0.1e1 / t667;
t688 = t652 * t657;
t687 = t662 * t669;
t686 = t664 * t671;
t685 = t666 * t673;
t684 = t668 * t669;
t683 = t670 * t671;
t682 = t672 * t673;
t674 = rSges(2,2) * g(3);
t681 = (((rSges(2,2) * t640 - t702) * t669 + (rSges(2,1) * t640 + t674) * t663) * m(2) + ((-t640 * rSges(3,3) + t703 * g(3)) * t669 + (-t703 * t640 - t701) * t663) * m(3)) * t650 / t668 ^ 2;
t680 = (((rSges(2,2) * t641 - t702) * t671 + (rSges(2,1) * t641 + t674) * t665) * m(2) + ((-t641 * rSges(3,3) + t704 * g(3)) * t671 + (-t704 * t641 - t701) * t665) * m(3)) * t651 / t670 ^ 2;
t679 = (((rSges(2,2) * t642 - t702) * t673 + (rSges(2,1) * t642 + t674) * t667) * m(2) + ((-t642 * rSges(3,3) + t705 * g(3)) * t673 + (-t705 * t642 - t701) * t667) * m(3)) * t652 / t672 ^ 2;
t675 = 0.1e1 / pkin(2);
t1 = [-m(4) * g(1) + (-(t645 * t666 + t648 * t682) * t688 - (t644 * t664 + t647 * t683) * t689 - (t643 * t662 + t646 * t684) * t690) * t700 + ((-t645 * t685 - t648 * t672) * t679 + (-t644 * t686 - t647 * t670) * t680 + (-t643 * t687 - t646 * t668) * t681 + (t643 * t693 + t644 * t692 + t645 * t691) * m(3)) * t675; -m(4) * g(2) + (-(t645 * t682 - t648 * t666) * t688 - (t644 * t683 - t647 * t664) * t689 - (t643 * t684 - t646 * t662) * t690) * t700 + ((-t672 * t645 + t648 * t685) * t679 + (-t670 * t644 + t647 * t686) * t680 + (-t668 * t643 + t646 * t687) * t681 + (-t646 * t693 - t647 * t692 - t648 * t691) * m(3)) * t675; (-m(4) - (3 * t649)) * g(3);];
taugX  = t1;
