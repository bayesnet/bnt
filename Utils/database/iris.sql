/*
 Navicat Premium Data Transfer

 Source Server         : growlithe
 Source Server Type    : MySQL
 Source Server Version : 80012
 Source Host           : localhost:3306
 Source Schema         : uci_database

 Target Server Type    : MySQL
 Target Server Version : 80012
 File Encoding         : 65001

 Date: 16/08/2019 10:52:54
*/

SET NAMES utf8mb4;
SET FOREIGN_KEY_CHECKS = 0;

-- ----------------------------
-- Table structure for iris
-- ----------------------------
DROP TABLE IF EXISTS `iris`;
CREATE TABLE `iris` (
  `id` int(10) NOT NULL AUTO_INCREMENT COMMENT 'id 主键',
  `sepal_length` decimal(10,2) DEFAULT NULL COMMENT '萼片长度',
  `sepal_width` decimal(10,2) DEFAULT NULL COMMENT '萼片宽度',
  `petal_length` decimal(10,2) DEFAULT NULL COMMENT '花瓣长度',
  `petal_width` decimal(10,2) DEFAULT NULL COMMENT '花瓣宽度',
  `class_name` varchar(20) CHARACTER SET utf8 COLLATE utf8_bin DEFAULT NULL COMMENT '类别名称',
  `status` int(1) DEFAULT NULL COMMENT '数据状态 1 有效 0 失效',
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=151 DEFAULT CHARSET=utf8 COLLATE=utf8_bin;

-- ----------------------------
-- Records of iris
-- ----------------------------
BEGIN;
INSERT INTO `iris` VALUES (1, 5.10, 3.50, 1.40, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (2, 4.90, 3.00, 1.40, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (3, 4.70, 3.20, 1.30, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (4, 4.60, 3.10, 1.50, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (5, 5.00, 3.60, 1.40, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (6, 5.40, 3.90, 1.70, 0.40, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (7, 4.60, 3.40, 1.40, 0.30, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (8, 5.00, 3.40, 1.50, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (9, 4.40, 2.90, 1.40, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (10, 4.90, 3.10, 1.50, 0.10, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (11, 5.40, 3.70, 1.50, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (12, 4.80, 3.40, 1.60, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (13, 4.80, 3.00, 1.40, 0.10, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (14, 4.30, 3.00, 1.10, 0.10, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (15, 5.80, 4.00, 1.20, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (16, 5.70, 4.40, 1.50, 0.40, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (17, 5.40, 3.90, 1.30, 0.40, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (18, 5.10, 3.50, 1.40, 0.30, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (19, 5.70, 3.80, 1.70, 0.30, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (20, 5.10, 3.80, 1.50, 0.30, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (21, 5.40, 3.40, 1.70, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (22, 5.10, 3.70, 1.50, 0.40, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (23, 4.60, 3.60, 1.00, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (24, 5.10, 3.30, 1.70, 0.50, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (25, 4.80, 3.40, 1.90, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (26, 5.00, 3.00, 1.60, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (27, 5.00, 3.40, 1.60, 0.40, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (28, 5.20, 3.50, 1.50, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (29, 5.20, 3.40, 1.40, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (30, 4.70, 3.20, 1.60, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (31, 4.80, 3.10, 1.60, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (32, 5.40, 3.40, 1.50, 0.40, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (33, 5.20, 4.10, 1.50, 0.10, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (34, 5.50, 4.20, 1.40, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (35, 4.90, 3.10, 1.50, 0.10, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (36, 5.00, 3.20, 1.20, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (37, 5.50, 3.50, 1.30, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (38, 4.90, 3.10, 1.50, 0.10, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (39, 4.40, 3.00, 1.30, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (40, 5.10, 3.40, 1.50, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (41, 5.00, 3.50, 1.30, 0.30, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (42, 4.50, 2.30, 1.30, 0.30, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (43, 4.40, 3.20, 1.30, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (44, 5.00, 3.50, 1.60, 0.60, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (45, 5.10, 3.80, 1.90, 0.40, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (46, 4.80, 3.00, 1.40, 0.30, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (47, 5.10, 3.80, 1.60, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (48, 4.60, 3.20, 1.40, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (49, 5.30, 3.70, 1.50, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (50, 5.00, 3.30, 1.40, 0.20, 'Iris-setosa', 1);
INSERT INTO `iris` VALUES (51, 7.00, 3.20, 4.70, 1.40, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (52, 6.40, 3.20, 4.50, 1.50, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (53, 6.90, 3.10, 4.90, 1.50, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (54, 5.50, 2.30, 4.00, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (55, 6.50, 2.80, 4.60, 1.50, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (56, 5.70, 2.80, 4.50, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (57, 6.30, 3.30, 4.70, 1.60, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (58, 4.90, 2.40, 3.30, 1.00, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (59, 6.60, 2.90, 4.60, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (60, 5.20, 2.70, 3.90, 1.40, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (61, 5.00, 2.00, 3.50, 1.00, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (62, 5.90, 3.00, 4.20, 1.50, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (63, 6.00, 2.20, 4.00, 1.00, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (64, 6.10, 2.90, 4.70, 1.40, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (65, 5.60, 2.90, 3.60, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (66, 6.70, 3.10, 4.40, 1.40, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (67, 5.60, 3.00, 4.50, 1.50, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (68, 5.80, 2.70, 4.10, 1.00, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (69, 6.20, 2.20, 4.50, 1.50, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (70, 5.60, 2.50, 3.90, 1.10, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (71, 5.90, 3.20, 4.80, 1.80, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (72, 6.10, 2.80, 4.00, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (73, 6.30, 2.50, 4.90, 1.50, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (74, 6.10, 2.80, 4.70, 1.20, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (75, 6.40, 2.90, 4.30, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (76, 6.60, 3.00, 4.40, 1.40, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (77, 6.80, 2.80, 4.80, 1.40, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (78, 6.70, 3.00, 5.00, 1.70, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (79, 6.00, 2.90, 4.50, 1.50, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (80, 5.70, 2.60, 3.50, 1.00, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (81, 5.50, 2.40, 3.80, 1.10, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (82, 5.50, 2.40, 3.70, 1.00, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (83, 5.80, 2.70, 3.90, 1.20, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (84, 6.00, 2.70, 5.10, 1.60, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (85, 5.40, 3.00, 4.50, 1.50, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (86, 6.00, 3.40, 4.50, 1.60, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (87, 6.70, 3.10, 4.70, 1.50, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (88, 6.30, 2.30, 4.40, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (89, 5.60, 3.00, 4.10, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (90, 5.50, 2.50, 4.00, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (91, 5.50, 2.60, 4.40, 1.20, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (92, 6.10, 3.00, 4.60, 1.40, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (93, 5.80, 2.60, 4.00, 1.20, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (94, 5.00, 2.30, 3.30, 1.00, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (95, 5.60, 2.70, 4.20, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (96, 5.70, 3.00, 4.20, 1.20, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (97, 5.70, 2.90, 4.20, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (98, 6.20, 2.90, 4.30, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (99, 5.10, 2.50, 3.00, 1.10, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (100, 5.70, 2.80, 4.10, 1.30, 'Iris-versicolor', 1);
INSERT INTO `iris` VALUES (101, 6.30, 3.30, 6.00, 2.50, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (102, 5.80, 2.70, 5.10, 1.90, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (103, 7.10, 3.00, 5.90, 2.10, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (104, 6.30, 2.90, 5.60, 1.80, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (105, 6.50, 3.00, 5.80, 2.20, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (106, 7.60, 3.00, 6.60, 2.10, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (107, 4.90, 2.50, 4.50, 1.70, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (108, 7.30, 2.90, 6.30, 1.80, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (109, 6.70, 2.50, 5.80, 1.80, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (110, 7.20, 3.60, 6.10, 2.50, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (111, 6.50, 3.20, 5.10, 2.00, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (112, 6.40, 2.70, 5.30, 1.90, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (113, 6.80, 3.00, 5.50, 2.10, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (114, 5.70, 2.50, 5.00, 2.00, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (115, 5.80, 2.80, 5.10, 2.40, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (116, 6.40, 3.20, 5.30, 2.30, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (117, 6.50, 3.00, 5.50, 1.80, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (118, 7.70, 3.80, 6.70, 2.20, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (119, 7.70, 2.60, 6.90, 2.30, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (120, 6.00, 2.20, 5.00, 1.50, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (121, 6.90, 3.20, 5.70, 2.30, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (122, 5.60, 2.80, 4.90, 2.00, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (123, 7.70, 2.80, 6.70, 2.00, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (124, 6.30, 2.70, 4.90, 1.80, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (125, 6.70, 3.30, 5.70, 2.10, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (126, 7.20, 3.20, 6.00, 1.80, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (127, 6.20, 2.80, 4.80, 1.80, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (128, 6.10, 3.00, 4.90, 1.80, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (129, 6.40, 2.80, 5.60, 2.10, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (130, 7.20, 3.00, 5.80, 1.60, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (131, 7.40, 2.80, 6.10, 1.90, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (132, 7.90, 3.80, 6.40, 2.00, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (133, 6.40, 2.80, 5.60, 2.20, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (134, 6.30, 2.80, 5.10, 1.50, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (135, 6.10, 2.60, 5.60, 1.40, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (136, 7.70, 3.00, 6.10, 2.30, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (137, 6.30, 3.40, 5.60, 2.40, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (138, 6.40, 3.10, 5.50, 1.80, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (139, 6.00, 3.00, 4.80, 1.80, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (140, 6.90, 3.10, 5.40, 2.10, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (141, 6.70, 3.10, 5.60, 2.40, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (142, 6.90, 3.10, 5.10, 2.30, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (143, 5.80, 2.70, 5.10, 1.90, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (144, 6.80, 3.20, 5.90, 2.30, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (145, 6.70, 3.30, 5.70, 2.50, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (146, 6.70, 3.00, 5.20, 2.30, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (147, 6.30, 2.50, 5.00, 1.90, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (148, 6.50, 3.00, 5.20, 2.00, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (149, 6.20, 3.40, 5.40, 2.30, 'Iris-virginica', 1);
INSERT INTO `iris` VALUES (150, 5.90, 3.00, 5.10, 1.80, 'Iris-virginica', 1);
COMMIT;

SET FOREIGN_KEY_CHECKS = 1;
