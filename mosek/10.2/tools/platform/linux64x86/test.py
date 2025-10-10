import mosek

# 创建一个Mosek任务环境
with mosek.Task() as task:
    # 简单地设置一些参数（这里只是示例）
    numvar = 2
    task.appendvars(numvar)
    # 打印任务信息
    print(task.getnumvar())

